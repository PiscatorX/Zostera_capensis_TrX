#! /usr/bin/env python
from Bio import SeqIO
import collections
import  sqlite3
import argparse
import pprint
import time
import sys
import os





class DB_Connect(object):
    
    def __init__(self, db_name, contigs= None, fix_id =  True, outfile = None):

        """

          Load contigs into SQLite database

        """
        
        self.db_name = db_name
        self.conx = sqlite3.connect(self.db_name)

        self.contigs  = contigs
        self.fix_id  = fix_id
        self.outfile = outfile if outfile else contigs.replace(".fasta","") +  "_fixed.fasta"  
                            

        
    def get_contig(self): 
        

        
        
        if not self.contigs:
            raise argparse.ArgumentTypeError('provide fasta contig file. See -h/--help')
        
        if self.fix_id:
            self.tsv_fp =  open('.'.join([self.outfile , "tsv"]), 'w') 
            self.contigs_newfname_fp  = open(self.outfile, 'w')         
            sys.stderr.write("{}\n{}\n".format(self.tsv_fp.name, self.contigs_newfname_fp.name))
            
        self.init_contigTable()
        contig_records = SeqIO.parse(self.contigs, "fasta")
        id_counts =  collections.defaultdict(int)
        contig_data = {}
        for contig in contig_records:
            id_counts[contig.id]=+1
            contig_id = contig.id + str(id_counts[contig.id]) if id_counts[contig.id] != 1 else contig.id
            descr = contig.description.split(" ",1)[1:][0]
            if self.fix_id:
                contig.id = contig_id
                contig.description = ''
                contig_data = dict([field.split("=") for field in descr.split(' ',1)]) 
                SeqIO.write(contig, self.contigs_newfname_fp,  "fasta")
                self.contigs_newfname_fp.flush()
                print("{}\t{}".format(contig.id, descr),file=self.tsv_fp, flush = True) 
            contig_data.update({'id': contig_id, 'length': str(len(contig))})
            cols = ', '.join(contig_data.keys())
            values = ', '.join([ "'"+val+"'" for val in contig_data.values() ])
            sql = "INSERT INTO contigs({}) VALUES ({})".format(cols, values )
            try:
                self.conx.execute(sql)
            except sqlite3.IntegrityError:
                break

        self.conx.commit()
        if self.fix_id:
            self.contigs_newfname_fp.close()
            self.tsv_fp.close()


            
            
    def init_contigTable(self):
        
        contigs_table = """
        CREATE TABLE IF NOT EXISTS contigs(
        id VARCHAR(20) PRIMARY KEY,
        flag INT,   
        len  INT, 
        path VARCHAR,
        length INT)"""
        self.conx.execute(contigs_table)

        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="""Connect to SQlite database""")
    parser.add_argument('contigs', help ='fasta file with contigs')
    parser.add_argument('-d','--db_name', help ='SQLite database filename', required = True)    
    parser.add_argument('-f','--fix-id', dest = "fix_id", action="store_true", help ='remove description from fasta def line')
    parser.add_argument('-o','--outfile', help = "Output filename for fixed contig and tsv filename", required = True)
    args = parser.parse_args()   
    db_connect = DB_Connect(args.db_name, args.contigs, outfile=args.outfile)
    db_connect.get_contig() 
