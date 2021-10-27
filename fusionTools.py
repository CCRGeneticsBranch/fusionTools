#!/usr/bin/env python
import re
import os,subprocess
import pandas as pd
import sys
import argparse
import tempfile
import threading
import json
from classes import *
from gtfparse import read_gtf
from Bio import SeqIO
from Bio.Seq import Seq
from pyfaidx import Fasta
from datetime import datetime
from dataclasses import dataclass
       
class FusionClassifier(threading.Thread):
    
    def __init__(self, genome, fusion_list, fusion_cancer_genes, fusion_cancer_pairs, cancer_genes, pfam_file):
        threading.Thread.__init__(self)
        self._genome = genome
        self._fusion_list = fusion_list        
        self._fusion_cancer_genes = fusion_cancer_genes
        self._fusion_cancer_pairs = fusion_cancer_pairs
        self._cancer_genes = cancer_genes
        self._pfam_file = pfam_file        
        self._results = {}
        self._detailed_results = {}
    
    
    def run(self):
        
        for key in self._fusion_list:
            row = key.split(":")
            left_symbol = row[0]
            right_symbol = row[1]            
            left_chr = row[2]
            right_chr = row[4]
            left_position = int(row[3])
            right_position = int(row[5])
            sample = row[6]
            left_fusion_cancer_gene = left_symbol in self._fusion_cancer_genes
            left_cancer_gene = left_symbol in self._cancer_genes
            right_fusion_cancer_gene = right_symbol in self._fusion_cancer_genes
            right_cancer_gene = right_symbol in self._cancer_genes
            
            left_gene = Gene(self._genome, left_symbol, left_chr, left_fusion_cancer_gene, left_cancer_gene)
            cancer_pair = (left_symbol + " " + right_symbol) in self._fusion_cancer_pairs
            right_gene = Gene(self._genome, right_symbol, right_chr, right_fusion_cancer_gene, right_cancer_gene)
            
            event = FusionEvent(self._genome, left_gene, left_position, right_gene, right_position, cancer_pair, self._pfam_file)
            (rep_result, fuse_peps, left_results, right_results, left_trans_info, right_trans_info) = event.process()
            gene_info = [left_fusion_cancer_gene, right_fusion_cancer_gene, left_cancer_gene, right_cancer_gene]
            self._results[key] = [self._fusion_list[key], rep_result, gene_info, fuse_peps, left_results, right_results, left_trans_info, right_trans_info]
    
    def getResult(self):
        return self._results
        
        
def main(args):
    #User input:    
    gtf_file = args.gtf.strip()
    fasta_file = args.fasta.strip()
    canonical_trans_file= args.canonical_trans_file.strip()
    in_file = args.input.strip()
    out_file = args.output.strip()
    pfam_file = args.pfam_file.strip()
    domain_file = args.domain_file.strip()
    num_threads = args.threads
    fusion_cancer_file = args.fusion_cancer_gene_list.strip()
    cancer_gene_file = args.cancer_gene_list.strip()
    isoform_expression_file = None
    if args.isoform_expression_file != None:
        isoform_expression_file = args.isoform_expression_file.strip()
        logging.info("Expression:" + isoform_expression_file)
    gene_bed_file = gtf_file.replace("gtf.gz", "genes.bed")
    logger = logging.getLogger()    
    if not os.path.exists(gene_bed_file):
        logger.info("Gene bed file not exists, generating one")
        gtf = read_gtf(gtf_file)
        gene_gtf = gtf[gtf['feature']=="gene"]
        gene_gtf = gene_gtf[["seqname","start","end","strand","gene_id","gene_name","level","gene_status"]]
        gene_gtf.to_csv(gene_bed_file, sep ='\t', index=None)    
    logging.info("GTF:" + gtf_file)
    logging.info("FASTA:" + fasta_file)
    logging.info("Input:" + in_file)    
    start=datetime.now()
    #prepare GTF, canonical list and cancer gene list
    genome = Genome(gtf_file, gene_bed_file, fasta_file, canonical_trans_file, domain_file, isoform_expression_file)
    threads = []
    avail_threads = os.cpu_count()
    logging.info("threads assigned:" + str(num_threads))
    logging.info("total cpus:" + str(avail_threads))
    if avail_threads < num_threads:
        num_threads = avail_threads - 1
    logging.info("total " + str(num_threads) + " are used")
    in_list = pd.read_csv(in_file, delimiter = "\t")
    fusion_list = {}
    #combine fusion callers
    for index, row in in_list.iterrows():
        if 1==2:
           target_gene = "FOXO1"
           if row[0] != target_gene and row[1] != target_gene:
              continue
        s = ":"
        rc = "NA"
        if "SpanReadCount" in row:
            rc = row["SpanReadCount"]
        key = s.join(map(str,row[0:7]))
        value = {row["Tool"]:rc}
        if key in fusion_list:
            fusion_list[key].append({row["Tool"]:rc})
        else:
            fusion_list[key] = [{row["Tool"]:rc}]
    chunks = []
    total_events = len(fusion_list)
    num_in_chunk = (int)(total_events/(num_threads-1))
    logger.info("total events:" + str(total_events))
    logger.info("num_in_chunk:" + str(num_in_chunk))
    i = 0
    chunk = {}
    #split data into chunks for multithreading
    for key in fusion_list:
        value = fusion_list[key]
        chunk[key] = value
        i = i + 1
        if i > num_in_chunk:
            chunks.append(chunk)
            chunk = {}
            i = 0
    chunks.append(chunk)
    fusion_cancer_genes = {}
    fusion_cancer_pairs = {}
    cancer_genes = {}
    pfam_file = pfam_file
    with open(fusion_cancer_file) as f:
        for line in f:
            (lg, rg) = line.split("\t")
            lg = lg.strip()
            rg = rg.strip()
            fusion_cancer_genes[lg] = ''
            fusion_cancer_genes[rg] = ''
            fusion_cancer_pairs[lg + " " + rg] = ''
    f.close
    with open(cancer_gene_file) as f:
        for line in f:
            cancer_genes[line.strip()] = ''
    f.close
    init_time = datetime.now()-start
    start=datetime.now()
    fusionClassifiers = []
    for c in chunks:
        fusionClassifier = FusionClassifier(genome, c, fusion_cancer_genes, fusion_cancer_pairs, cancer_genes, pfam_file)
        fusionClassifier.start()
        fusionClassifiers.append(fusionClassifier)
    
    for f in fusionClassifiers:
        f.join()
    
    #output results. We use tab seperated text
    sep = "\t"
    of = open(out_file,"w")
    header = ["left_gene", "right_gene", "left_chr", "right_chr", "left_position", "right_position", "sample_id", "tools", "type", "tier", \
    "left_region", "right_region", "left_trans", "right_trans", "left_fusion_cancer_gene", "right_fusion_cancer_gene", "left_cancer_gene", "right_cancer_gene", "fusion_proteins", "left_trans_info", "right_trans_info"]
    of.write(sep.join(header) + "\n")
    for f in fusionClassifiers:
        results = f.getResult()
        for key in results:   
            key_str = key.replace(":", sep)
            tools, rep_result, gene_info, fuse_peps, left_results, right_results, left_trans_info, right_trans_info = results[key]
            rep_str = sep.join(map(str, [rep_result["type"],  rep_result["tier"], rep_result["left_location"], rep_result["right_location"], rep_result["left_trans"], rep_result["right_trans"]]))
            gene_info_str = sep.join(map(lambda x: "Y" if x else "N", gene_info))
            of.writelines(key_str + "\t" + json.dumps(tools) + "\t" + rep_str + "\t" + gene_info_str + "\t" + json.dumps(fuse_peps) + "\t" + json.dumps(left_trans_info) + "\t" + json.dumps(right_trans_info) + "\n")
    of.close()
    process_time = datetime.now()-start    
    logger.info(init_time)
    logger.info(process_time)
    
    
    
script_dir = os.path.dirname(os.path.abspath(__file__))
avail_threads = os.cpu_count() -1
parser = argparse.ArgumentParser(description='Classify fusion types.')
parser.add_argument("--gtf", "-g", metavar="GTF file", default=script_dir + "/data/gencode.v38lift37.annotation.sorted.gtf.gz", help="GTF file")
parser.add_argument("--fasta", "-f", metavar="Genome FASTA file", default=script_dir + "/data/hg19.fasta.gz", help="Genome FASTA file")
parser.add_argument("--isoform_expression_file", "-m", metavar="Isoform expression file in RSEM format", help="Isoform expression file in RSEM format")
parser.add_argument("--canonical_trans_file", "-n", metavar="Canonical transcript list", default=script_dir + "data/canonical_transcripts_with_ensembl.txt", help="Conoical transcript list (default: %(default)s)")
parser.add_argument("--fusion_cancer_gene_list", "-u", metavar="Fusion cancer gene pair list", default=script_dir + "data/sanger_mitelman_pairs.txt", help="Fusion cancer gene pair list (default: %(default)s)")
parser.add_argument("--cancer_gene_list", "-c", metavar="Cancer gene list", default=script_dir + "data/clinomics_gene_list.txt", help=" (default: %(default)s)Cancer gene list")
parser.add_argument("--pfam_file", "-p", metavar="Pfam domain file", default=script_dir + "/PfamDB", help="Pfam domain file (default: %(default)s)")
parser.add_argument("--input", "-i", metavar="Fusion file", required=True, help="[Fusion input file]")
parser.add_argument("--output", "-o", metavar="output file", help="[output file]", required=True)
parser.add_argument("--domain_file", "-d", metavar="Pfam domain file", default=script_dir + "/data/gencode.v38lift37.domains.tsv", help="[Pfam domain file (default: %(default)s)]")
parser.add_argument("--threads", "-t", metavar="(Number of threads)", type=int, default=avail_threads, help="[Number of threads (default: %(default)s)]")
try:
    args = parser.parse_args()
except:
    parser.print_help()
    sys.exit(0)

main(args)