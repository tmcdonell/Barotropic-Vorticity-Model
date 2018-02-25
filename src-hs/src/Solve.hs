{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE RebindableSyntax    #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE ViewPatterns        #-}
-- |
-- Module      : Solve
-- Copyright   : [2018] Trevor L. McDonell
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <tmcdonell@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--

module Solve
  where

import Type

import Data.Array.Accelerate
import Data.Array.Accelerate.Control.Lens


-- Laplacian operator
--
-- src-c/main.cpp:122
--
laplacian
    :: forall a. Fractional a
    => Exp a
    -> Exp a
    -> Acc (Field a)
    -> Acc (Field a)
laplacian dx dy = stencil update boundary
  where
    update :: Stencil3x3 a -> Exp a
    update ((_,t,_)
           ,(l,c,r)
           ,(_,b,_)) = (l + r + t + b - 4.0 * c) / (dx * dy)

    -- out-of-bounds elements treated as zero
    boundary :: Boundary (Field a)
    boundary = function (const 0)


-- Implements the basic time-stepping of the vorticity equation:
--
-- dz/dt = - ( d psi dx * dzdy - d psi dy * dzdx )
--
-- src-c/main.cpp:140
--
advance
    :: forall a. Fractional a
    => Exp a
    -> Exp a
    -> Exp a                -- time step
    -> Exp a                -- gradient of the corriolis parameter
    -> Acc (Field a)        -- stream function psi at time n
    -> Acc (Field a)        -- vorticity at time (n-1)
    -> Acc (Field a)        -- vorticity at time n
    -> Acc (Field a)        -- vorticity at time (n+1)
advance dx dy dt beta psi prev curr =
  zipWith3 (\u v w -> u - dt * (v + w)) prev (arakawa dx dy psi curr) (stencil s clamp psi)
  where
    s :: Stencil3x3 a -> Exp a
    s (_,(l,_,r),_) = beta * (r-l) / (2*dx)


-- Arakawa Jacobian
--
-- src-c/main.cpp:187
--
arakawa
    :: forall a. Fractional a
    => Exp a
    -> Exp a
    -> Acc (Field a)    -- first variable (stream function)
    -> Acc (Field a)    -- second variable (vorticity)
    -> Acc (Field a)
arakawa dx dy var1 var2 = stencil2 update wrap var1 wrap var2
  where
    update :: Stencil3x3 a -> Stencil3x3 a -> Exp a
    update v1 v2 =
      let ((a1,b1,c1),
           (d1,_ ,f1),
           (g1,h1,i1)) = v1

          ((a2,b2,c2),
           (d2,_ ,f2),
           (g2,h2,i2)) = v2
          --
          j1  = (f1-d1) * (h2-b2) - (f2-d2) * (h1-b1)
          j2  = f1*(i2-c2) - d1*(g2-a2) - h1*(i2-g2) - b1*(c2-a2)
          j3  = h2*(i1-g1) - b2*(c1-a1) - f2*(i1-c1) - d1*(g1-a1)
      in
      (1/3) * (j1+j2+j3) / (4.0*dx*dy)


-- Solve the Poisson equation using successive over-relaxation (SOR)
-- The equation is of the form:
--
--   div . grad(var) = rhs
--
-- src-c/main.cpp:264
--
poisson
    :: forall a. (Ord a, Fractional a, FromIntegral Int a)
    => Exp a
    -> Exp a
    -> Exp a
    -> Exp a
    -> Acc (Field a)
    -> Acc (Field a)
    -> Acc (Field a)
poisson relTol omega dx dy var0 rhs
  = view _3
  $ awhile
      (view _2)                                     -- while this returns True...
      (\st ->                                       -- ...keep executing this update function
         let var       = st^._3       -- old approximation
             var'      = step var     -- updated value
             (i, cont) = continue (st^._1) var var'
         in
         lift (i, cont, var'))
      (lift (unit 0, unit (constant True), var0))   -- ...beginning with this initial state
  where
    Z :. height :. width  = unlift (shape var0)

    lIMIT = width * height
    dx2   = dx*dx
    dy2   = dy*dy
    tau   = 0.5 * dx2*dx2 / (dx2 + dy2)

    residual :: Acc (Field a) -> Acc (Field a) -> Acc (Scalar a)
    residual var var'
      = sum
      $ flatten
      $ zipWith (\v v' -> abs (1 - v/v')) var var'

    continue :: Acc (Scalar Int) -> Acc (Field a) -> Acc (Field a) -> (Acc (Scalar Int), Acc (Scalar Bool))
    continue i var var'
      = unzip
      $ zipWith (\n r -> let n' = n+1
                             r' = r / fromIntegral (width*height)
                         in
                         lift (n', n' < lIMIT && r' > relTol))
                i
                (residual var var')

    step :: Acc (Field a) -> Acc (Field a)
    step var = stencil2 update boundary var clamp rhs
      where
        update :: Stencil3x3 a -> Stencil3x3 a -> Exp a
        update v1 v2 =
          let ((_,t,_),
               (l,c,r),
               (_,b,_))     = v1
              (_,(_,x,_),_) = v2    -- just the current value of the RHS
          in
          c + omega * tau / dx2 * (r + l - 2*c)
            + omega * tau / dy2 * (t + b - 2*c)
            - omega * tau * x

        -- Periodic in the x-direction, zero at the y-boundary
        -- The index given to us is out-of-bounds it an least one dimension.
        --
        boundary :: Boundary (Field a)
        boundary = function $ \(unlift -> Z :. y :. x) ->
          if y < 0 || y >= height
            then 0
            else if x < 0
                   then var ! index2 y (width + x)
                   else var ! index2 y (x - width)

