-- Flubud calculates the Rohlin distance between alphabetical strings,
-- according to a certain procedure of partitioning.
--
-- For details, see:
-- Burioni R , Scalco R , Casartelli M (2011) Rohlin Distance and the Evolution of
-- Influenza A Virus: Weak Attractors and Precursors. PLoS ONE 6(12): e27924.
-- doi:10.1371/journal.pone.0027924
--
-- Copyright (C) 2012  Riccardo Scalco
-- 
-- This program is free software: you can redistribute it and/or modify
-- it under the terms of the GNU General Public License as published by
-- the Free Software Foundation, either version 3 of the License, or
-- (at your option) any later version.
-- 
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
-- GNU General Public License for more details.
-- 
-- You should have received a copy of the GNU General Public License
-- along with this program.  If not, see <http://www.gnu.org/licenses/>.
--
-- Email: riccardo.scalco@gmail.com


module Flubud.Core
( getSequenceList
, getPartition
, getPositionsAA
, getAtoms
, spanPair
, matRohlin
, distRohlin
, muFlat
, entropy
, minComRef
, sortFirst
, toEdges
, fromEdges
, reduction
, distRohlinR
, matRohlinR
) where


import Data.List ((\\), intersect, sort, concat)
import qualified Data.Map as Map (fromListWith, Map, elems)
import qualified Data.Graph as Graph
import qualified Data.Set as Set
import qualified Data.Tree as Tree


type Part = [[Int]]


-- | Return the up-diagonal matrix of all pair distances between the partitions on the
-- list @p:l@ according the measure @mu@.
-- > let p = [[[1,2],[3,4]],[[1],[2,3,4]],[[1],[2],[3,4]]]
-- > matRohlin muFlat p == [[3.295836866004329,1.3862943611198906],[1.9095425048844386]]
matRohlinR           :: ([Int] -> Double) -> [Part] -> [[Double]]
matRohlinR mu [p]    = []
matRohlinR mu (p:l)  = (map (distRohlinR mu p) l): matRohlinR mu l

-- | Return the Rohlin distance between partitions @x@ and @y@.
-- > let a = [[1,3,4],[2,5]]
-- > let b = [[1],[2,3,4],[5]]
-- > distRohlin muFlat a b == 5.2053793708887675
distRohlinR         :: ([Int] -> Double) -> Part -> Part -> Double
distRohlinR mu x y  = 2 * entropy mu (minComRef xR yR) - entropy mu xR - entropy mu yR
   where (xR,yR)    = reduction (x,y)



-- |
-- 
-- reduction :: (Part,Part) -> (Part,Part)
-- reduction ([],[]) = ([],[])
reduction (a,b)  =
    let a_base   = toEdges [concat a]
        b_base   = toEdges [concat b]
        amb      = Set.difference a_set b_set
        bma      = Set.difference b_set a_set
        a_set    = toEdges a
        b_set    = toEdges b
        bounds   = (1,length $ concat a)
        a_red    = fromEdges bounds $ Set.difference a_base bma
        b_red    = fromEdges bounds $ Set.difference b_base amb
    in (a_red, b_red)




-- | Return the set of edges from the partition.
-- > toEdges [[1,2],[3],[4,5],[6]] == Set.fromList [(1,2),(4,5)]
-- > toEdges [[1,2],[3,4,5],[6]]   == Set.fromList [(1,2),(3,4),(4,5)]
toEdges     :: Part -> Set.Set (Int,Int)
toEdges xs  = Set.fromList $ foldl (\acc x -> acc ++ zip (init x) (tail x)) [] xs




-- | Return the partition from the set of edges and the bounds of the
-- partition.
-- > (fromEdges (1,6) $ toEdges [[1,2],[3],[4,5],[6]]) == [[1,2],[3],[4,5],[6]] 
fromEdges :: (Int,Int) -> Set.Set (Int,Int) -> Part 
fromEdges bounds s =
    let forest = Graph.components graph
        graph = Graph.buildG bounds (Set.toList s)
    in sortPart [Tree.flatten i | i <- forest] 



-- > let s = [">first","hor","se",">second","ele","phant"]
-- > getSequenceList s == ["horse","elephant"]
-- > let s = [">first",">second","hor","se",">third","ele","phant"]
-- > getSequenceList s == ["","horse","elephant"]
getSequenceList     :: [[Char]] -> [[Char]]
getSequenceList []  = []
getSequenceList xs  =
                    let (ys,zs) = spanPair (\x y -> head y /= '>') xs
                    in concat (tail ys) : getSequenceList zs



-- | Sort the partition.
-- > sortPart [[2,6],[5,4],[3,1]] == [[1,3],[2,6],[4,5]]
sortPart :: Part -> Part
sortPart x = sort [sort i | i <- x]



-- | Return the partition associated to the sequence @l@. The parameter @s@
-- represent the accepted interval between AA.
-- > getPartition 0 "aababbac" == [[1],[2],[3],[4],[5],[6],[7],[8]]
-- > getPartition 1 "aababbac" == [[1,2],[3],[4],[5,6],[7],[8]]
-- > getPartition 2 "aababbac" == [[1,2,4],[3,5,6],[7],[8]]
-- > getPartition 3 "aababbac" == [[1,2,4,7],[3,5,6],[8]]
getPartition       :: Int -> [Char] -> Part
getPartition s []  = []
getPartition s l   = sort $ concat [getAtoms s (sort a) | a <- Map.elems (getPositionsAA l)]



-- | Return a Map where key are the AA in the sequence and values their positions.
-- > getPositionsAA "aababbac"
--   fromList [('a',[7,4,2,1]),('b',[6,5,3]),('c',[8])]
getPositionsAA    :: [Char] -> Map.Map Char [Int]
getPositionsAA l  = Map.fromListWith (++) $ map (\(k,v) -> (k,[v])) (zip l [1..length l])



-- | Return the atoms from the list of all the position labels. It is assumed that
-- the labels are ordered. The parameter @s@ represent the accepted difference
-- between labels.
-- > getAtoms 0 [1,2,4,5,8] == [[1],[2],[4],[5],[8]]
-- > getAtoms 1 [1,2,4,5,8] == [[1,2],[4,5],[8]]
-- > getAtoms 2 [1,2,4,5,8] == [[1,2,4,5],[8]]
-- > getAtoms 3 [1,2,4,5,8] == [[1,2,4,5,8]]
getAtoms        :: Int -> [Int] -> Part
getAtoms s []   = []
getAtoms s [x]  = [[x]]
getAtoms s xs   =
                let (ys,zs) = spanPair (\x y -> y-x <= s) xs
                in ys : getAtoms s zs



-- | 'spanPair', applied to a predicate @p@ and a list @xs@, returns a tuple where
-- first element is longest prefix of @xs@ of elements that
-- satisfy @p@ with their precedent, and second element is the remainder of the list.
-- > spanPair (\x y -> abs(x-y) <= 1) [1,2,4,6,7]  == ([1,2],[4,6,7])
-- > spanPair (\x y -> abs(x-y) <= 2) [1,2,4,8,13] == ([1,2,4],[8,13])
-- > spanPair (\x y -> x+y<=6) [1,5,1,5,13]        == ([1,5,1,5],[13])
-- > spanPair (\x y -> x+y<=6) [1,5,1,8,13]        == ([1,5,1],[8,13])
-- > spanPair (\x y -> x+y<=6) [1,5,1,5,10,1,2]    == ([1,5,1,5],[10,1,2])
spanPair                    :: (a -> a -> Bool) -> [a] -> ([a],[a])
spanPair _ xs@[]            =  (xs, xs)
spanPair _ xs@[x]           =  (xs,[])
spanPair p xs@(x:xs')
         | p x (head xs')   =  let (ys,zs) = spanPair p xs' in (x:ys,zs)
         | otherwise        =  ([x],xs')



-- | Return the up-diagonal matrix of all pair distances between the partitions on the
-- list @p:l@ according the measure @mu@.
-- > let p = [[[1,2],[3,4]],[[1],[2,3,4]],[[1],[2],[3,4]]]
-- > matRohlin muFlat p == [[3.295836866004329,1.3862943611198906],[1.9095425048844386]]
matRohlin           :: ([Int] -> Double) -> [Part] -> [[Double]]
matRohlin mu [p]    = []
matRohlin mu (p:l)  = (map (distRohlin mu p) l): matRohlin mu l



-- | Return the Rohlin distance between partitions @x@ and @y@.
-- > let a = [[1,3,4],[2,5]]
-- > let b = [[1],[2,3,4],[5]]
-- > distRohlin muFlat a b == 5.2053793708887675
distRohlin         :: ([Int] -> Double) -> Part -> Part -> Double
distRohlin mu x y  = 2 * entropy mu (minComRef x y) - entropy mu x - entropy mu y



-- | Return the length of the atom @x@ as measure of @x@.
-- > muFlat [1,5,6] == 3.0
muFlat    :: [Int] -> Double
muFlat x  = fromIntegral (length x)



-- | Return the entropy of the partition @p@ according the measure @mu@.
-- > entropy muFlat [[1,2],[3]] == -1.3862943611198906
entropy       :: ([Int] -> Double) -> Part -> Double
entropy mu p  =
              let ent x = - mu x * log (mu x)
              in foldl (\acc x -> acc + ent x) 0 p


-- | Return the minimum common refinement between two partitions,
-- under the assumption that partitions are ordered.
-- > minComRef [[1,2],[3]] [[1],[2,3]] == [[1],[2],[3]]
-- > minComRef [[1,3,4],[2,5]] [[1],[2,3,4],[5]] == [[1],[2],[3,4],[5]]
minComRef             :: Part -> Part -> Part
minComRef [] []       = []
minComRef ([]:xs) ys  = minComRef xs ys
minComRef xs ([]:ys)  = minComRef xs ys
-- minComRef (x:xs) (x:ys) = x: minComRef xs ys
minComRef (x:xs) (y:ys)
                      | x == y            = x : minComRef xs ys
                      | head x == head y  = intersection : minComRef (compl_x:xs) (compl_y:ys)
                      | otherwise         = minComRef (sortFirst $ x:xs) (sortFirst $ y:ys)
                      where intersection  = intersect x y
                            compl_x       = x \\ intersection
                            compl_y       = y \\ intersection



-- | Sort a list under the assumption that only the first argument is not in
-- order.
-- > sortFirst [2,1,4,5,6] == [1,2,4,5,6]
-- > sortFirst [7,1,8,5] == [1,7,8,5] (the assumption does not hold here!)
-- > sortFirst [[2],[1,3],[4,5]] == [[1,3],[2],[4,5]]
sortFirst        :: (Ord a) => [a] -> [a]
sortFirst []     = []
sortFirst [x]    = [x]
sortFirst (x:y:xs)
                 | x > y     = y: sortFirst (x:xs)
                 | otherwise = x:y:xs


