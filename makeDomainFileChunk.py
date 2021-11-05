#!/usr/bin/python
import re,gzip
import os,subprocess
import pandas as pd
import json
import pandas as pd
from classes import *
        
def main(in_file, out_file, pfam):
    df = pd.read_csv(in_file, delimiter = "\t", header=None)
    with open(out_file, "w+") as of:
        for index, row in df.iterrows():
            sequence = row[1]
            domains = predictPfam(sequence, pfam)
            #domains = []
            of.write(row[0] + "\t" + sequence + "\t" + str(row[2]) + "\t" + json.dumps(domains) + "\n")
    of.close()

script_dir = os.path.dirname(os.path.abspath(__file__))
parser = argparse.ArgumentParser(description='Classify fusion types.')
parser.add_argument("--in_file", "-i", metavar="input protein FASTA file", help="Protein FASTA file")
parser.add_argument("--out_file", "-o", metavar="output text file", help="Output domain file")
parser.add_argument("--pfam", "-p", metavar="Pfam DB directory", default=script_dir + "/PfamDB", help="Output domain file")
args = parser.parse_args()
in_file = args.in_file.strip()
out_file = args.out_file.strip()
pfam = args.pfam.strip()
main(in_file, out_file, pfam)