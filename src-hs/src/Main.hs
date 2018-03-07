{-# LANGUAGE BangPatterns      #-}
{-# LANGUAGE CPP               #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
-- |
-- Module      : Main
-- Copyright   : [2018] Trevor L. McDonell
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <tmcdonell@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--

module Main where

import Type
import Solve

import Text.Printf
import Control.Lens
import System.IO
import Data.Text.Lazy                                               ( Text )
import Data.Text.Lazy.Read                                          as T
import qualified Data.Text.Lazy                                     as T
import qualified Data.Text.Lazy.IO                                  as T
import Prelude                                                      as P

import Data.Array.Accelerate                                        as A ( Z(..), (:.)(..) )
import qualified Data.Array.Accelerate                              as A
#ifdef ACCELERATE_LLVM_NATIVE_BACKEND
import Data.Array.Accelerate.LLVM.Native                            as CPU
#endif
#ifdef ACCELERATE_LLVM_PTX_BACKEND
import Data.Array.Accelerate.LLVM.PTX                               as PTX
#endif


input :: FilePath
input = "../data/barotropic_streamfunction_input.dat"

output :: Int -> FilePath
output n = printf "out/streamfunction_%04d.dat" n


main :: IO ()
main = do
  -- Read in grid data. Requires data to be whitespace separated in the same
  -- shape as the output array (i.e. number of lines -> number of rows in the
  -- matrix, number of fields per row -> number of columns).
  --
  streamFunction <- parse <$> T.readFile input

  let -- Simulation parameters and other crap
      --
      deltaX              = 200.0e3                           :: R
      deltaY              = 200.0e3                           :: R
      Z :. sizeY :. sizeX = A.arrayShape streamFunction
      beta                = 1.0e-11                           :: R
      relTol              = 1.0e-7                            :: R
      deltaT              = 15 * 60                           :: R
      forecastTime        = 144 * 3600                        :: R
      numTimeSteps        = P.round (forecastTime / deltaT)   :: Int
      omega               = 1.5                               :: R
      --
      dx                  = A.constant deltaX
      dy                  = A.constant deltaY
      dt                  = A.constant deltaT

  printf "input file name:         %s\n" input
  printf "num x grid points:       %d\n" sizeX
  printf "num y grid points:       %d\n" sizeY
  printf "time step:               %f seconds\n" deltaT
  printf "delta x:                 %f\n" deltaX
  printf "delta y:                 %f\n" deltaY
  printf "number of time steps:    %d\n" numTimeSteps
  printf "using omega value:       %f\n" omega
  printf "\n"

  let -- Initialise vorticity field
      --
      zeta0   = A.fill (A.shape psi0) 0
      psi0    = A.use streamFunction

      -- Take first time step
      --
      zeta1   = advance dx dy (2*dt) (A.constant beta) psi0 zeta0 zeta0
      psi1    = poisson (A.constant relTol) (A.constant omega) dx dy psi0 zeta0

      start   = CPU.run $ A.lift (zeta0, zeta1, psi1)

      -- Solve for the new stream function
      --
      !go     = CPU.runN (step dx dy (2*dt) (A.constant beta) (A.constant omega) (A.constant relTol))

      loop !i !st
        | i > numTimeSteps  = return ()
        | otherwise         = do
            let !st'  = go st
                path  = output i
            --
            printf "\rprogress: %3.0f %%" (100 * fromIntegral i / fromIntegral numTimeSteps :: Double)
            if i `rem` 12 == 0
              then do printf "\rwriting data: %s\n" path
                      dump path (st'^._3)
              else return ()
            --
            loop (i+1) st'

  -- main loop
  loop 0 start


-- For really large files we may want to read in the data multiple times to get
-- the width & height, so that we don't keep the entire file contents resident
-- in memory (UTF16 encoded).
--
parse :: Text -> Field R
parse txt =
  let
      width   = length (T.words (head (T.lines txt)))
      height  = length (T.lines txt)
      --
      ok (Right (v, rest))
        | T.null rest     = realToFrac v
        | otherwise       = error $ printf "parse error: unhandled input: %s" (T.unpack rest)
      ok (Left err)       = error $ printf "parse error: %s" err
  in
  A.fromList (Z :. height :. width) [ ok (signed double v) | v <- T.words txt ]

dump :: FilePath -> Field R -> IO ()
dump path arr =
  withFile path WriteMode $ \h ->
  let
      Z :. _ :. width = A.arrayShape arr

      go !_ []     = return ()
      go !i (x:xs) = do
        hPutStr h (printf "%g" x)
        if i < width-1
          then hPutChar h ' '  >> go (i+1) xs
          else hPutChar h '\n' >> go 0     xs
  in
  go 0 (A.toList arr)

