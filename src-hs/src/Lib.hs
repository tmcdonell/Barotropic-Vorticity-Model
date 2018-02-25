{-# LANGUAGE RebindableSyntax #-}
-- |
-- Module      : Lib
-- Copyright   : [2018] Trevor L. McDonell
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <tmcdonell@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--

module Lib (

  dotp

) where

import Data.Array.Accelerate

-- | A simple vector inner product
--
dotp :: Acc (Vector Double) -> Acc (Vector Double) -> Acc (Scalar Double)
dotp xs ys = fold (+) 0 ( zipWith (*) xs ys)
