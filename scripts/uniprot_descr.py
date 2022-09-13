#!/usr/bin/env  python
import itertools as it
import argparse


def getx(def_fileObj, tags):
    
    columns = ["ACC"]+tags
    for line in def_fileObj:
        
        idx, descx = line.strip().split("=",1)
        idx = [ idx2 for idx1 in  idx.split(" ",1) for idx2 in idx1.rsplit(" ",1) ]
        descx = [ descx2 for descx1 in descx.split("=") for descx2 in descx1.rsplit(" ",1)  ]
        descx.insert(0, idx.pop(-1))
        descx_iter = iter(descx)
        entry = dict(zip(descx_iter, descx_iter))
        entry.update(dict(zip(["ACC","ID","name"],[ x_i for x in idx for x_i in x.split("|") ][1:])))
        
        print("\t".join([entry[col] for col in columns ]))


        
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("def_file", type=argparse.FileType('r'))
    parser.add_argument('-t','--tags', nargs='+', required = True)
    args = parser.parse_args()
    getx(args.def_file, args.tags)    







        
    

        
