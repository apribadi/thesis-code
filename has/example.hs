module Main where

import Data.Array
import Data.List
import qualified Data.Foldable as Foldable
import qualified System.Random.MWC as Rand



-- a_ij, b_i, c_j
type Param = (Array (Int, Int) Double, Array Int Double, Array Int Double)

-- Lexicographic order, i.e. integer keys correspond to their binary
-- representation.
type Dist = Array Int Double

-- number of hidden and visible nodes
kk = 2
nn = 3

if' b x y = if b then x else y


parametrization :: Param -> Dist
parametrization (w, b, c) = array (0, 2^nn-1) [(i, psi (vs ! i) / z) | i <- [0, 2^nn-1]]
  where
    energy v h = 
        sum [if' (v ! i && h ! j) (w ! (i, j)) 0 | i <- [1..nn], j <- [1..kk]]
      + sum [if' (v ! i) (b ! i) 0 | i <- [1..nn]]
      + sum [if' (h ! j) (c ! j) 0 | j <- [1..kk]]
    psi v = sum [energy v h | h <- Foldable.toList hs]
    vs = bitvecs nn
    hs = bitvecs kk
    z = sum [energy v h | v <- Foldable.toList vs, h <- Foldable.toList hs]

-- least significant bit first
toBinary n = unfoldr (\k -> if (k == 0) then Nothing else Just (mod k 2, div k 2)) n
toBools n = map (== 1) (toBinary n)


-- nondet choice from lists
-- lexicographic order
choices :: [[a]] -> [[a]]
choices []     = [[]]
choices (x:xs) = [y:zs | y <- x, zs <- choices xs]

bitvecs :: Int -> Array Int (Array Int Bool)
bitvecs n = listArray (0, 2^n-1) [listArray (1, n) xs | xs <- choices $ replicate n [False, True]]


hausdorff xs ys f = max a b
  where
    a = maximum [minimum [f x y | y <- ys] | x <- xs]
    b = maximum [minimum [f x y | x <- xs] | y <- ys]

    
main = do
    --putStrLn $ show $ hausdorff [1,2,3] [3,4,5] (\x y -> abs (x - y))
    print $ choices $ replicate 5 [0,1]
