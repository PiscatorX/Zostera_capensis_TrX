from Bio import SwissProt
import pprint
import argparse


def getattr_swissprot(sprot_fileObj,  tags):
    
    i  = 1
    for record in SwissProt.parse(sprot_fileObj):
        print(getattr(record, 'gene_name')) 
        # attr_data = {}
        # for tag in tags:
        #     attr = getattr(record, tag)
        #     attr_data[tag] = attr
        #     if attr:
        #         attr_data[tag] = attr.split(";",1)[0].split("=")[1]
    
        # print(attr_data)
        i = i + 1
        if i == 50:
            break
        
        
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("sprot_dat", type=argparse.FileType('r'))
    parser.add_argument('-t','--tags', nargs='+', required = True)
    args = parser.parse_args()
    getattr_swissprot(args.sprot_dat, args.tags)    
