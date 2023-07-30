#!/usr/bin/env python
from Bio import AlignIO
import argparse
import pprint
import sys



def pfam_GF(pfam_seed, columns):
    
    records = AlignIO.parse(open(pfam_seed, encoding="ISO-8859-1"), "stockholm")
    for rec in records:
        print('\t'.join([ rec.general_features[k][0] for k  in columns ]))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "convert msa file")
    parser.add_argument('pfam_seed')
    parser.add_argument('-c','--columns', nargs='+', type = str, required = True)
    args = parser.parse_args()
    pfam_GF(args.pfam_seed, args.columns)
