#!/usr/bin/env  python
from Bio import SwissProt
import argparse
import pprint
import sys


def get_mappingx(swissprot_dat_file, accessionlist,  xref, column):

    sp_handle = open(swissprot_mapping_file)
    sp_db =  { acc:record for record in SwissProt.parse(sp_handle) for  acc in record.accessions }

    acc_list  = tuple([ line for line in accessionlist.read().strip().splitlines()])
    for acc in acc_list:
        sp_record = sp_db.get(acc.split("\t")[column-1], None)
        if sp_record:
            for entry in sp_record.cross_references:            
                if xref in entry: 
                    entry = list(entry)
                    entry.pop(entry.index(xref))
                    print("{}\t{}".format(acc, "\t".join(entry)))
        else:
            print (acc)

                    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("swissprot_dat_file")
    parser.add_argument('-m','--mappinglist', type=argparse.FileType('r'), required = True, help = "list file of accession numbers, one per line")
    parser.add_argument('-x','--xref', default = "GO")
    parser.add_argument('-c','--column', type=int, default = 2, help = "list file of accession numbers, one per line")
    args = parser.parse_args()
    get_mappingx(args.swissprot_dat_file, args.mappinglist, args.xref, args.column)    







        
    

        
