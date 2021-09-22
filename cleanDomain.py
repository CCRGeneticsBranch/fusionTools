#!/usr/bin/python
import re
import os,subprocess
import pandas as pd
import json
from classes import *
        
def main():
    #gtf_file = "/data/Clinomics/Ref/khanlab/GTF/ucsc.hg19_star.gtf"
    gtf_file = "/data/Clinomics/Ref/khanlab/GTF/gencode.v19.annotation.gtf"
    fasta_file = "/data/Clinomics/Ref/khanlab/ucsc.hg19.fasta"
    pfam = "/data/khanlab/projects/hsienchao/oncogenomics/fusion/fusionTools/PfamDB"
    genome = Genome(gtf_file, fasta_file, None)    
    df = pd.read_csv("ensembl_domain.old.tsv", delimiter = "\t")
    print("trans\tseq\tprotein_length\tdomains")
    for index, row in df.iterrows():
        t = row["TRANS"]
        aa_seq = str(row["AA_SEQ"])
        l = len(aa_seq)
        try:
            domains = json.loads(row["DOMAIN"].replace('""', '"'))
            nd = []        
            for domain in domains:
                s = int(domain["start_pos"].strip())
                e = int(domain["end_pos"].strip())
                n = domain["name"].strip()
                a = domain["hint"]["Accession"].strip()
                reg = re.search('>(.*)<', a, re.IGNORECASE)
                if reg:
                    a = reg.group(1)
                d = domain["hint"]["Description"].strip()
                nd.append([s,e,n,a,d])
            print(t + "\t" + str(aa_seq) + "\t" + str(l) + "\t" + json.dumps(nd))
        except:
            tx = Transcript(genome, genome, t)
            seq = tx.getCodingSequence()
            domains = predictPfam(seq, pfam)
            print(t + "\t" + seq + "\t" + str(l) + "\t" + json.dumps(domains))
    
main()