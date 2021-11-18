from gtfparse import read_gtf
import sys, argparse

parser = argparse.ArgumentParser(description='Make gene region bed file form GTF.')
parser.add_argument("--gtf", "-i", metavar="GTF file")
parser.add_argument("--bed", "-o", metavar="BED file")

args = parser.parse_args()

gtf_file = args.gtf.strip()
bed_file = args.bed.strip()

gtf = read_gtf(gtf_file)
gene_gtf = gtf[gtf['feature']=="gene"]
gene_gtf = gene_gtf[["seqname","start","end","strand","gene_id","gene_name"]]
gene_gtf.to_csv(bed_file, sep ='\t', index=None, header=False)
