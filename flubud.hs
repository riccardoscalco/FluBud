{-# LANGUAGE TemplateHaskell #-}

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

import Data.List ((\\), intersect, sort, concat)
import qualified Data.Map as Map (fromListWith, Map, elems)
import Options
import System.IO (readFile)
import Flubud.Core

defineOptions "MainOptions" $ do
    intOption "optSpace" "s" 1
        "The accepted difference between labels (optional, default to 1)."
    stringOption "optFastaFile" "fasta" ""
        "The fasta file to read as input (required)."

main :: IO ()
main = runCommand $ \opts args -> do
    if null (optFastaFile opts)
        then error "The --fasta flag is required."
        else return ()
    contents <- readFile $ optFastaFile opts
    let ps   = [getPartition (optSpace opts) l | l <- getSequenceList $ lines contents]
        mat  = matRohlinR muFlat ps
    mapM_ print mat
