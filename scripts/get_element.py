#!/usr/bin/env python
import argparse
import sys

def cutx(myfile, column, separator):
    
    for line in myfile.read().splitlines(): 
        elements = line.split(separator)
        n = len(elements)
        if n >= column:
            print(elements[column - 1])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "finiky splitter")
    parser.add_argument('myfile', type=argparse.FileType('r'))
    parser.add_argument('-c','--column', type = int, default = 2)
    parser.add_argument('-s','--separator', default = "\t")
    args = parser.parse_args()
    cutx(args.myfile, args.column, args.separator)

    
