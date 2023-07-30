#!/usr/bin/env  python
import itertools as it
import argparse
import sys


def getx(def_fileObj, separator, out_separator):

    for line in def_fileObj:
        
        sep_cols = line.strip().split(separator)
        sep_cols = [ col_j for col_i in  sep_cols for col_j in col_i.split(separator) ]
        print(out_separator.join(sep_cols))
        

#/pipexsplit.py  data.enterprise/Blast_top_hits/EvigeneX.blastx |cut -f 1,3  > data.enterprise/Blast_top_hits/EvigeneX.blastx.ids        
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("def_file", type=argparse.FileType('r'))
    parser.add_argument('-s','--separator', default = "|")
    parser.add_argument('-o','--out_separator', default = "\t")
    args = parser.parse_args()
    getx(args.def_file, args.separator, args.out_separator)    







        
    

        
