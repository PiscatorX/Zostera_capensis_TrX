import pandas as pd
import sys
import csv

infname = sys.argv[1]
data = csv.DictReader(open(infname),  delimiter = "\t")
csv_writer = csv.DictWriter(open(infname+".rst", "w"), fieldnames = ['USA','Length'], delimiter = "\t")
csv_writer.writeheader()
i = 0
for row in data:
    USA = row['USA']    
    row['USA']  = USA.split(":")[2]
    _  = [ row.pop(k, None) for k in  ['Database', 'Name', 'Accession', 'Type', '%GC', 'Organism', 'Description' ]]
    csv_writer.writerow(row)
