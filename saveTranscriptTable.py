from gtfparse import read_gtf
import pandas as pd
import sys, argparse

def remove_ensembl_version(s):
    a = s.split(".", 1)
    return a[0]
    
parser = argparse.ArgumentParser(description='Make gene region bed file form GTF.')
parser.add_argument("--gtf", "-g", metavar="GTF file")
parser.add_argument("--out", "-o", metavar="Output file")

args = parser.parse_args()

gtf_file = args.gtf.strip()
out_file = args.out.strip()

gtf = read_gtf(gtf_file)
trans_gtf = gtf.loc[gtf['feature']=="transcript"]
trans_gtf["transcript_id"] = trans_gtf["transcript_id"].apply(lambda x: remove_ensembl_version(x))
trans_gtf.to_csv(out_file, sep="\t", index=False)