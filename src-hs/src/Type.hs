-- |
-- Module      : Type
-- Copyright   : [2018] Trevor L. McDonell
-- License     : BSD3
--
-- Maintainer  : Trevor L. McDonell <tmcdonell@cse.unsw.edu.au>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--

module Type
  where

import Data.Array.Accelerate
import Data.Array.Accelerate.Linear.V2


type R       = Double
type Field a = Array DIM2 a

-- Notes regarding the reference implementation (modelGrid.cpp)
--
--  - Implicitly returns 0 when read out-of-bounds
--  - Also keeps the delta_x and delta_y values
--
type Grid a = Field (V2 a)

