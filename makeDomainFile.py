#!/usr/bin/python
import re,gzip
import os,subprocess
import pandas as pd
import json
from classes import *
        
def main(in_file, out_file, pfam):
    with gzip.open(in_file, "rt") as handle:
        fasta_sequences = SeqIO.parse(handle,'fasta')
        with open(out_file, "w+") as of:
            #of.write("trans\tseq\tprotein_length\tdomains\n")
            for fasta in fasta_sequences:
                name, sequence = fasta.id, str(fasta.seq)
                id_search = re.search(".*\|(ENST.*)\|ENSG.*", name)
                if id_search:
                    trans_id = id_search.group(1);
                    trans_id = remove_ensembl_version(trans_id)
                    length = len(sequence)
                    #domains = predictPfam(sequence, pfam)
                    domains = []
                    of.write(trans_id + "\t" + sequence + "\t" + str(length) + "\t" + json.dumps(domains) + "\n")
        of.close()

script_dir = os.path.dirname(os.path.abspath(__file__))
parser = argparse.ArgumentParser(description='Classify fusion types.')
parser.add_argument("--in_file", "-i", metavar="input protein FASTA file", default=script_dir + "/data/gencode.v36lift37.pc_translations.fa.gz", help="Protein FASTA file")
parser.add_argument("--out_file", "-o", metavar="output text file", default=script_dir + "/data/gencode.v36lift37.pc_translations.tsv", help="Output domain file")
parser.add_argument("--pfam", "-p", metavar="Pfam DB directory", default=script_dir + "/PfamDB", help="Output domain file")
args = parser.parse_args()
in_file = args.in_file.strip()
out_file = args.out_file.strip()
pfam = args.pfam.strip()
main(in_file, out_file, pfam)