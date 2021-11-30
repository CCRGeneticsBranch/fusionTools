#!/usr/bin/python
import re
import os,subprocess,io
import warnings
from Bio import BiopythonWarning
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter('ignore', BiopythonWarning)
import pandas as pd
import yaml
import sys
import argparse
import tempfile
import threading
import json
import pysam
import logging
from gtfparse import read_gtf
from Bio import SeqIO
from Bio.Seq import Seq
from pyfaidx import Fasta
from datetime import datetime
from dataclasses import dataclass
from pybedtools import BedTool

STATUS_CODING = "STATUS_CODING"
STATUS_INTERGENIC = "STATUS_INTERGENIC"
STATUS_UTR = "STATUS_UTR"
STATUS_NO_CDS = "STATUS_NO_CDS"
STATUS_NO_EXONS = "STATUS_NO_EXONS"
STATUS_INTRON_SPLICE_SITE_PROXIMAL = "STATUS_INTRON_SPLICE_SITE_PROXIMAL"
STATUS_INTRON_SPLICE_SITE_DISTAL = "STATUS_INTRON_SPLICE_SITE_DISTAL"
TYPE_IN_FRAME = "in-frame"
TYPE_OUT_OF_FRAME = "out-of-frame"
TYPE_RIGHT_INTACT = "right gene intact"
TYPE_NO_PROTIEN = "no protein"
TYPE_NO_ANNOTATION = "no annotation"

def predictPfam(fuse_pep, pfam_file):
    tmp_in = os.path.basename(tempfile.NamedTemporaryFile().name)
    out_data = ">fuse\n" + fuse_pep
    f = open(tmp_in, "w")        
    f.write(out_data)
    f.flush()
    f.close()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    #os.system("hmmscan --domtblout " + tmp_out + " -E 1e-5 --domE 1e-5 " + self._pfam_file + " " + tmp_in)
    cmd = "perl " + script_dir + "/PfamScan/pfam_scan.pl -fasta " + tmp_in + " -dir " + pfam_file
    stream = os.popen(cmd)
    domains = []
    for line in stream.readlines():
        line = line.strip()
        if line[0:1] != "#":
            #print(line)
            tokens = line.split()            
            if len(tokens) > 5:
                domains.append([int(tokens[1]), int(tokens[2]), tokens[6], tokens[5], tokens[6]])
    f.close()
    os.remove(tmp_in)
    return domains

def remove_ensembl_version(s):
    a = s.split(".", 1)
    return a[0]
    
class Genome:
    def __init__(self, gtf_file, gene_bed_file, fasta_file, canonical_trans_file=None, domain_file=None, isoform_expression_file=None):
        self._gtf = pysam.TabixFile(gtf_file)
        self._gene_bedtool_obj = BedTool(gene_bed_file)
        self._gene_bed = pd.read_csv(gene_bed_file, delimiter = "\t")
        self._gene_bed.columns = ['seqname','start','end','strand','gene_id','gene_name']
        gene_ids = self._gene_bed["gene_id"].apply(lambda x: remove_ensembl_version(x))
        self._gene_bed["gene_id"] = gene_ids
        self._gene_bed["gene_name"] = self._gene_bed["gene_name"].str.upper()
        self._gene_bed.set_index("gene_name")

        #self._gtf = read_gtf(gtf_file)
        self._sequences = Fasta(fasta_file, one_based_attributes=True)
        if canonical_trans_file != None:
            self._canonical_trans = pd.read_csv(canonical_trans_file, delimiter = "\t")
            self._canonical_trans["Symbol"] = self._canonical_trans["Symbol"].str.upper()
            self._canonical_trans.set_index("Symbol")
        if domain_file != None:
            self._domains = pd.read_csv(domain_file, delimiter = "\t")
            self._domains.set_index("trans")
        self._isoform_expression = None
        if isoform_expression_file != None:
            self._isoform_expression = pd.read_csv(isoform_expression_file, delimiter = "\t")
            isoform_ids = self._isoform_expression["transcript_id"].apply(lambda x: remove_ensembl_version(x))
            self._isoform_expression["transcript_id"] = isoform_ids
            self._isoform_expression.set_index("transcript_id")
            #print(self._isoform_expression.head(5))
        #gene_ids = self._gtf["gene_id"].apply(lambda x: remove_ensembl_version(x))
        #trans_ids = self._gtf["transcript_id"].apply(lambda x: remove_ensembl_version(x))
        #self._gtf["gene_id"] = gene_ids
        #self._gtf["transcript_id"] = trans_ids
        #self._gtf["gene_name"] = self._gtf["gene_name"].str.upper()
        #self._gtf.set_index("gene_name")
    
    @property
    def gene_bedtool_obj(self):
        return self._gene_bedtool_obj
        
    @property
    def gene_bed(self):
        return self._gene_bed
        
    @property
    def gtf(self):
        return self._gtf
    @property
    def domains(self):
        return self._domains
    @property
    def sequences(self):
        return self._sequences
        
    @property
    def canonical_trans(self):
        return self._canonical_trans
    
    @property
    def isoform_expression(self):
        return self._isoform_expression
        
    def getSequence(self, chromosome, start, end, oneBased=True):
        start = int(start)
        end = int(end)
        if oneBased:
            start = start - 1
        #print(chromosome + ":" + str(start) + "-" + str(end))
        return str(self._sequences[chromosome][start:end])

@dataclass
class Exon:
    start_pos: int
    end_pos: int
    exon_number: int
    type: str
    
    def __iter__(self):
        yield 'start_pos', self.start_pos
        yield 'end_pos', self.end_pos
        yield 'exon_number', self.exon_number
        yield 'type', self.type
        
    def __str__(self):
        return json.dumps(dataclasses.asdict(self))

@dataclass
class FusionElement:
    transcript: None
    status: str
    location: str
    exon_number: str
    sequence: str
    extend_info: None

    
class Transcript:

    def __init__(self, genome, gene, transcript_id):
        self._genome = genome
        self._gene = gene
        self._transcript_id = transcript_id
        self._expression = "NA"
        self._is_canonical = (transcript_id == gene.canonical_transcript_id);
        self._is_mane = (transcript_id == gene.mane_transcript_id);
        self._trans_gtf = gene.gtf[gene.gtf["transcript_id"]==transcript_id]
        if genome.isoform_expression is not None:
            exp = genome.isoform_expression[genome.isoform_expression["transcript_id"]==transcript_id]            
            if not exp.empty:
                self._expression = exp.iloc[0]["TPM"]
        #self._trans_gtf = gene.gtf[gene.gtf["transcript_id"].str.contains(transcript_id)]
        if not self._trans_gtf.empty:
            self._start = min(self._trans_gtf["start"])
            self._end = max(self._trans_gtf["end"])
            self._exon_info = []            
            self._domains = "[]"
            self._protein_length = 0
            self._strand = self._trans_gtf.iloc[0]["strand"]
            if hasattr(genome, 'domains'):
                tid = transcript_id.split(".", 1)[0]
                domains_df = genome.domains[genome.domains["trans"]==tid]
                if not domains_df.empty:
                    self._domains = domains_df.iloc[0]["domains"]
                    self._protein_length = int(domains_df.iloc[0]["protein_length"]) + 1
            self._trans_gtf = self._trans_gtf[self._trans_gtf['feature'] != "gene"]
            #print(self._strand)
    
    def getCodingSequence(self, do_translate=True):
        if self._trans_gtf.empty:
            return ""
        cds = self.getFeature()
        tx_seq = ""
        strand = "+"
        for ei, exon in cds.iterrows():
            seq = exon["seqname"]
            start = exon["start"]
            end = exon["end"]
            strand = exon["strand"]            
            cds_seq = self._genome.getSequence(seq, start, end)            
            tx_seq = tx_seq + cds_seq
        if strand == "-":
            tx_seq = str(Seq(tx_seq).reverse_complement())
        if do_translate:
            tx_seq = Seq(tx_seq).translate()
        return str(tx_seq)
   
    def getFusionSequence(self, breakpoint, upstream=True, intron_max=100):        
        location = ""
        tx_seq = ""
        exon_number = "NA"
        extend_info = {"offset":0, "seq":""}
        extend_cutoff = 8
        if self._trans_gtf.empty:
            return FusionElement(self, STATUS_NO_EXONS, location, exon_number, "", extend_info)
        status = STATUS_CODING
        # breakpoint not in gene region
        if breakpoint < self._start:
            location = "upstream" if self._strand=="+" else "downstream"
            status = STATUS_INTERGENIC
        if breakpoint > self._end:
            location = "downstream" if self._strand=="+" else "upstream"
            status = STATUS_INTERGENIC        
        exons = self.getFeature('exon')
        if exons.empty:
            return FusionElement(self, STATUS_NO_EXONS, location, exon_number, "", extend_info)
        exons = exons.astype({'exon_number': 'int32'})
        max_en = max(exons["exon_number"])
        cdses = self.getFeature('CDS')
        cds_start = 0
        cds_end = 0
        if not cdses.empty:
            cds_start = min(cdses["start"])
            cds_end = max(cdses["end"])
        previous_start = -1
        previous_end = -1
        previous_en = -1        
        #self._exon_info = []
        #loop exons to identify breakpoint is in exon/intron and exon number
        for ei, exon in exons.iterrows():
            start = exon["start"]
            end = exon["end"]
            en = exon["exon_number"]
            #gencode GTF has correct exon number. For refseq, use the following to reverse the exon number
            #if self._strand == "-":
            #    en = int(max_en) - int(en) + 1
            #exon info
            if not cdses.empty:
                #if CDS
                if start >= cds_start and end <= cds_end:
                    self._exon_info.append(dict(Exon(start, end, en, "CDS")))
                #if start/stop codon in exon
                if cds_start > start and cds_start < end:                    
                    self._exon_info.append(dict(Exon(start, cds_start - 1, en, "UTR5" if self._strand == "+" else "UTR3")))
                    self._exon_info.append(dict(Exon(cds_start, end, en, "CDS")))
                if cds_end > start and cds_end < end:
                    self._exon_info.append(dict(Exon(start, cds_end, en, "CDS")))
                    self._exon_info.append(dict(Exon(cds_end + 1, end, en, "UTR3" if self._strand == "+" else "UTR5")))
                #if whole exon is UTR
                if cds_start > start and cds_start > end:
                    self._exon_info.append(dict(Exon(start, end, en, "UTR5" if self._strand == "+" else "UTR3")))
                if cds_end < start and cds_end < end:
                    self._exon_info.append(dict(Exon(start, end, en, "UTR3" if self._strand == "+" else "UTR5")))
            else:
                self._exon_info.append(dict(Exon(start, end, en, "exon")))
            #print(str(en))
            # breakpoint in exon
            if breakpoint >= start and breakpoint <= end:
                location = "exon"
                exon_number = "exon" + str(en)
            # breakpoint in intron
            if breakpoint > previous_end and breakpoint < start:
                location = "intron"
                exon_number = "exon" + str(previous_en) + "-exon" + str(en) if en > previous_en else "exon" + str(en) + "-exon" + str(previous_en)
                if (breakpoint - previous_end) > intron_max and (start - breakpoint) > intron_max: 
                    status = STATUS_INTRON_SPLICE_SITE_DISTAL
            previous_start = start
            previous_end = end
            previous_en = en
        
        
        if cdses.empty:
            status = STATUS_NO_CDS        
        if status != STATUS_INTERGENIC: 
            if breakpoint < cds_start and status != STATUS_NO_CDS and exon_number != "NA":
                utr_type = "5UTR" if self._strand == "+" else "3UTR"
                location = utr_type
                status = STATUS_UTR            
            if breakpoint > cds_end and status != STATUS_NO_CDS and exon_number != "NA":
                utr_type = "3UTR" if self._strand == "+" else "5UTR"
                location = utr_type
                status = STATUS_UTR            
        
        #we don't need to calculate the fusion sequence
        if status == STATUS_INTRON_SPLICE_SITE_DISTAL or status == STATUS_NO_CDS:
            return FusionElement(self, status, location, exon_number, "", extend_info)
            
        #loop CDSs to get sequence
        previous_start = -1
        previous_end = -1
        #print("first part") if bool(upstream) == bool(self._strand == "+") else print("second part")
                
        for ei, cds in cdses.iterrows():
            seq = cds["seqname"]
            start = cds["start"]
            end = cds["end"]
            cds_seq = ""
            # we want the first part of the transcript
            if bool(upstream) == bool(self._strand == "+"):
                if breakpoint >= start:                    
                    if breakpoint > end:
                        #whole CDS included
                        cds_seq = self._genome.getSequence(seq, start, end)
                        tx_seq = tx_seq + cds_seq
                    else:
                        # breakpoint in CDS
                        cds_seq = self._genome.getSequence(seq, start, breakpoint)
                        location = "CDS"
                        tx_seq = tx_seq + cds_seq
                        if end - breakpoint <= extend_cutoff and end != breakpoint:
                            extend_info["offset"] = end - breakpoint
                            extend_info["seq"] = self._genome.getSequence(seq, breakpoint + 1, end)
                        break
                #check if in intron
                else:
                    if breakpoint > previous_end:
                        if breakpoint - previous_end > intron_max:
                            return FusionElement(self, STATUS_INTRON_SPLICE_SITE_DISTAL, location, exon_number, "", extend_info)
                        status = STATUS_INTRON_SPLICE_SITE_PROXIMAL
                        cds_seq = self._genome.getSequence(seq, previous_end+1, breakpoint)
                        extend_info["offset"] = -1*len(cds_seq)
                        extend_info["seq"] = cds_seq
                        tx_seq = tx_seq + cds_seq
                        break
            # we want the second part of the transcript
            else:
                if breakpoint <= end:
                    if breakpoint < start:                        
                        if breakpoint > previous_end and status != STATUS_UTR:
                            #in intron
                            status = STATUS_INTRON_SPLICE_SITE_DISTAL
                            if start - breakpoint < intron_max:
                                status = STATUS_INTRON_SPLICE_SITE_PROXIMAL
                                cds_seq = self._genome.getSequence(seq, breakpoint, end)
                                tx_seq = tx_seq + cds_seq
                                extend_info["offset"] = -1*(start - breakpoint)
                                extend_info["seq"] = cds_seq[:(start - breakpoint)]
                        else:
                            #whole CDS included
                            cds_seq = self._genome.getSequence(seq, start, end)
                            tx_seq = tx_seq + cds_seq
                    else:
                        # breakpoint in CDS
                        cds_seq = self._genome.getSequence(seq, breakpoint, end)
                        location = "CDS"
                        tx_seq = tx_seq + cds_seq
                        if breakpoint - start <= extend_cutoff and start != breakpoint:
                            extend_info["offset"] = breakpoint - start
                            extend_info["seq"] = self._genome.getSequence(seq, start, breakpoint - 1)
            previous_start = start
            previous_end = end
        if self._strand == "-":
            tx_seq = str(Seq(tx_seq).reverse_complement())
            if extend_info["offset"] != 0:
                extend_info["seq"] = str(Seq(extend_info["seq"]).reverse_complement())
        #print(extend_info)
        return FusionElement(self, status, location, exon_number, tx_seq, extend_info)
    
    @property
    def is_canonical(self):
        return self._is_canonical
    
    @property
    def is_mane(self):
        return self._is_mane
    
    def getFeature(self, feature="CDS"):
        if feature == "CDS":
            try:
                return self._trans_gtf.loc[(self._trans_gtf['feature'] == "CDS") | (self._trans_gtf['feature'] == "stop_codon")].sort_values(by=['start'])
            except Exception as inst:
                print("Exception in " + self._transcript_id + " " + self._gene.symbol)
        return self._trans_gtf.loc[self._trans_gtf['feature'] == feature].sort_values(by=['start'])
        #return self._trans_gtf.query('feature == "' + feature + '"').sort_values(by=['start'])    
    
    @property
    def exonInfo(self):
        return self._exon_info
    
    @property
    def strand(self):
        return self._strand
        
    @property
    def domains(self):
        return self._domains
    
    @property
    def expression(self):
        return self._expression
        
    @property
    def protein_length(self):
        return self._protein_length
        
    @property
    def transcript_id(self):
        return self._transcript_id
        
class Gene:
    
    def __init__(self, genome, symbol, chromosome, position, fusion_cancer_gene=False, cancer_gene=False):
        #Arriba might output multiple symbols
        self._symbol_label = symbol
        symbols = symbol.split(",")
        if len(symbols) > 1:
            min_dist = sys.maxsize
            current_symbol = symbol
            for s in symbols:
                search=re.search('(.*)\((.*)\)', s)
                if search:
                    sb = search.group(1)
                    dist = search.group(2)
                    if dist.isnumeric():
                        dist = int(dist)
                        if dist < min_dist:
                            current_symbol = sb
                            min_dist = dist
            symbol = current_symbol
        self._genome = genome        
        self._symbol = symbol
        self._chromosome = chromosome
        self._fusion_cancer_gene = fusion_cancer_gene
        self._cancer_gene = cancer_gene
        self._cancer_gene = cancer_gene
        self._canonical_transcript_id = ""
        self._mane_transcript_id = ""
        self._transcript_list = []
        g_df = genome.gene_bed.loc[(genome.gene_bed["gene_name"]==symbol.upper()) & (genome.gene_bed["seqname"]==chromosome)]
        #if cannot find symbol, try Ensembl accession
        if g_df.empty and symbol[0:4] == "ENSG":
            g_df = genome.gene_bed.loc[(genome.gene_bed["gene_id"]==symbol) & (genome.gene_bed["seqname"]==chromosome)]
        #if still empty, look for nearest gene
        if g_df.empty:
            new_symbol = self.getNearestSymbol(genome.gene_bedtool_obj, chromosome, position)
            if new_symbol != "":
                self._symbol_label = symbol + "(" + new_symbol + ")"
                symbol = new_symbol
                self._symbol = symbol
                g_df = genome.gene_bed.loc[(genome.gene_bed["gene_name"]==symbol.upper()) & (genome.gene_bed["seqname"]==chromosome)]
            else:
                return

        #tbx = pysam.TabixFile("data/gencode.v38lift37.annotation.sorted.gtf.gz")
        s = io.StringIO("\n".join(genome.gtf.fetch(g_df.iloc[0].seqname, g_df.iloc[0].start, g_df.iloc[0].end, multiple_iterators=True)))
        #s = io.StringIO("\n".join(tbx.fetch(g_df.iloc[0].seqname, g_df.iloc[0].start, g_df.iloc[0].end)))
        logging.getLogger().setLevel(logging.ERROR)
        self._gtf = read_gtf(s)
        self._gtf["gene_id"] = self._gtf["gene_id"].apply(lambda x: remove_ensembl_version(x))
        self._gtf["transcript_id"] = self._gtf["transcript_id"].apply(lambda x: remove_ensembl_version(x))
        self._gtf["gene_name"] = self._gtf["gene_name"].str.upper()
        self._gtf.set_index("gene_name")        
        self._gtf = self._gtf.loc[(self._gtf["gene_name"]==symbol.upper()) | (self._gtf["gene_id"]==symbol.upper())]
        #self._gtf = genome.gtf[genome.gtf["gene_name"]==symbol.upper()]
        #if self._gtf.empty and symbol[0:4] == "ENSG":
            #self._gtf = genome.gtf[genome.gtf["gene_id"]==symbol]
        
        can_txs = genome.canonical_trans[genome.canonical_trans["Symbol"]==symbol.upper()]        
        if not can_txs.empty:
            if can_txs.iloc[0]["Canonical"] == "Y":
                self._canonical_transcript_id = can_txs.iloc[0]["Ensembl"]
            if can_txs.iloc[0]["MANE"] == "Y":
                self._mane_transcript_id = can_txs.iloc[0]["Ensembl"]    
        self._transcript_list = list(set(self._gtf.loc[self._gtf["feature"]=="transcript"]["transcript_id"].tolist()))
    
    def getNearestSymbol(self, bedtools_obj, chromosome, position):
        bp_bed = BedTool(chromosome+"\t" + str(position) + "\t" + str(position),from_string=True)
        symbol = ""
        nearby = bp_bed.closest(bedtools_obj, d=True,stream=True)
        for gene in nearby:
            symbol = gene[-2]
        return symbol
    
    @property
    def symbol_label(self):
        return self._symbol_label
        
    @property
    def gtf(self):
        return self._gtf
        
    @property
    def symbol(self):
        return self._symbol
    
    @property
    def chromosome(self):
        return self._chromosome
        
    @property
    def fusion_cancer_gene(self):
        return self._fusion_cancer_gene
    
    @property
    def cancer_gene(self):
        return self._cancer_gene
        
    @property
    def canonical_transcript_id(self):
        return self._canonical_transcript_id
        
    @property
    def mane_transcript_id(self):
        return self._mane_transcript_id
    @property
    def transcript_list(self):
        return self._transcript_list
        

class FusionEvent:
    
    def __init__(self, genome, left_gene, left_position, right_gene, right_position, fusion_cancer_pairs, pfam_file):
        self._genome = genome
        self._left_gene = left_gene
        self._left_position = left_position
        self._right_gene = right_gene
        self._right_position = right_position
        self._fusion_cancer_pairs = fusion_cancer_pairs
        self._pfam_file = pfam_file
        s = ":"
        self._key = s.join(map(str,[left_gene.symbol,right_gene.symbol,left_gene.chromosome,left_position,right_gene.chromosome,right_position]))   
        self._rep_result = None
        self._fuse_peps = {}
        self._all_results = {}
        self._left_results = {}
        self._right_results = {}
        self._left_trans_info = {}
        self._right_trans_info = {}
    
    def __str__(self):
        sep = "\t"
        rep_str = sep.join(map(str, [self._left_gene.symbol, self._right_gene.symbol, self._left_position, self._right_position, \
        self._rep_result["type"],  self._rep_result["tier"], self._rep_result["left_trans"], self._rep_result["right_trans"]]))
        
        return rep_str + "\t" + json.dumps(self._fuse_peps)
        #return rep_str + "\t" + json.dumps(fuse_pep_dic) + "\t" + json.dumps(self._left_results) + "\t" + json.dumps(self._right_results)
    
    @property
    def left_gene(self):
        return self._left_gene
    
    @property
    def right_gene(self):
        return self._right_gene
        
    @property
    def left_position(self):
        return self._left_position
    
    @property
    def right_position(self):
        return self._right_position
    
    @property
    def fusion_cancer_pairs(self):
        return self._fusion_cancer_pairs
    
    @property
    def key(self):
        return self._key
        
    @property
    def rep_result(self):
        return self._rep_result
    
    @property
    def fuse_peps(self):
        return self._fuse_peps
        
    def getTier(self, type):
        if self._fusion_cancer_pairs:
            if type == TYPE_IN_FRAME:
                return 1.1
            else:
                return 1.2
        if self._left_gene.fusion_cancer_gene or self._right_gene.fusion_cancer_gene:
            if type == TYPE_IN_FRAME:
                return 2.1
            else:
                return 2.2
        if self._left_gene.cancer_gene or self._right_gene.cancer_gene:
            if type == TYPE_IN_FRAME or type == TYPE_RIGHT_INTACT:
                if self._left_gene.cancer_gene and self._right_gene.cancer_gene:
                    return 2.3
                else:
                    return 3.1
            return 3.2
        if type == TYPE_IN_FRAME:
            return 4.1
        if type == TYPE_RIGHT_INTACT:
            return 4.2
        return 4.3
        
    def process(self):
        self._all_results = {}
        self._left_results = {}
        self._right_results = {}
        self.left_info = {}
        self.right_info = {}
        # get transcript fusion part first
        for left_tx in self._left_gene.transcript_list:
            left_trans = Transcript(self._genome, self._left_gene, left_tx)
            #if left_trans.protein_length <=10:
            #    continue
            self._left_results[left_tx] = left_trans.getFusionSequence(self._left_position, True)
            self._left_trans_info[left_tx] = {"chr": self._left_gene.chromosome, "strand": left_trans.strand, "is_canonical": left_trans.is_canonical, "is_mane": left_trans.is_mane, "exon_info": left_trans.exonInfo, "expression": left_trans.expression, "seq_to_ss": self._left_results[left_tx].extend_info["seq"], "domains": left_trans.domains, "protein_length": left_trans.protein_length, "cdna":self._left_results[left_tx].sequence}
        for right_tx in self._right_gene.transcript_list:
            right_trans = Transcript(self._genome, self._right_gene, right_tx)
            #if right_trans.protein_length <=10:
            #    continue
            self._right_results[right_tx] = right_trans.getFusionSequence(self._right_position, False)
            self._right_trans_info[right_tx] = {"chr": self._right_gene.chromosome, "strand": right_trans.strand, "is_canonical": right_trans.is_canonical, "is_mane": right_trans.is_mane, "exon_info": right_trans.exonInfo,  "expression": right_trans.expression, "seq_to_ss": self._right_results[right_tx].extend_info["seq"], "domains": right_trans.domains, "protein_length": right_trans.protein_length, "cdna":self._right_results[right_tx].sequence}
        # get fused transcript for all combination
        fuse_peps = {}
        rep_tier = 99
        canonical_tier = 99
        # if no annotation found
        if len(self._left_gene.transcript_list) == 0 or len(self._right_gene.transcript_list) == 0:
            tier = self.getTier(TYPE_NO_ANNOTATION)
            self._rep_result = {"left_trans":"NA", "right_trans":"NA", "type": TYPE_NO_ANNOTATION, "tier": tier, \
                "fuse_nucl_position": 0, "fuse_pep_position" : 0, "left_location": "NA", "left_exon_number": "NA", "right_location": "NA", "right_exon_number": "NA"}
        #look for representative fusion: we use canonical transcripts. if there is no defined canonical transcripts, use the best tier pair
        fused_domains = {}
        for left_tx in self._left_trans_info.keys():
            is_left_canonical = (left_tx == self._left_gene.canonical_transcript_id)
            left_result = self._left_results[left_tx]
            for right_tx in self._right_trans_info.keys():
                key = left_tx + ":" + right_tx
                right_result = self._right_results[right_tx]
                is_right_canonical = (right_tx == self._right_gene.canonical_transcript_id)
                fuse_pep = ""
                type = TYPE_NO_PROTIEN
                fuse_nucl_seq = ""
                fuse_nucl_position = 0
                fuse_pep_position = 0
                #right gene intact. Ignore the left gene
                if (right_result.location == "upstream" or right_result.location == "5UTR") and right_result.sequence != "":
                    type = TYPE_RIGHT_INTACT
                    fuse_nucl_seq = right_result.sequence
                    fuse_pep = str(Seq(fuse_nucl_seq).translate())
                    fused_domains[fuse_pep] = self._right_trans_info[right_tx]["domains"]
                elif left_result.sequence != "" and right_result.sequence != "":
                    #if left_result.extend_info["offset"] != 0 and right_result.extend_info["offset"] != 0:
                    #    print(left_result.extend_info["seq"])
                    #    print(right_result.extend_info["seq"])                        
                    fuse_nucl_position = len(left_result.sequence) + 1
                    fuse_pep_position = int(fuse_nucl_position / 3) + 1
                    fuse_nucl_seq = left_result.sequence + right_result.sequence
                    fuse_pep = str(Seq(fuse_nucl_seq).translate())
                    first_stop = fuse_pep.find('*')
                    if first_stop == -1 or first_stop == len(fuse_pep) - 1:
                        type = TYPE_IN_FRAME
                    else:
                        type = TYPE_OUT_OF_FRAME
                        fuse_pep = fuse_pep[0:first_stop+1]
                tier = self.getTier(type)
                fusion_result = {"left_trans":left_tx, "right_trans":right_tx, "type": type, "tier": tier, \
                "fuse_nucl_position": fuse_nucl_position, "fuse_pep_position" : fuse_pep_position, \
                "left_location": left_result.location, "left_exon_number": left_result.exon_number, "right_location": right_result.location, "right_exon_number": right_result.exon_number}
                self._all_results[key] = fusion_result
                if fuse_pep != "":
                    if fuse_pep in fuse_peps:
                        fuse_peps[fuse_pep].append(fusion_result)
                    else:
                        fuse_peps[fuse_pep] = [fusion_result]
                if is_left_canonical and is_right_canonical:
                    canonical_tier = tier
                    self._rep_result = fusion_result
                if (not is_left_canonical) and (self._left_gene.canonical_transcript_id != ""):
                    continue
                if (not is_right_canonical) and (self._right_gene.canonical_transcript_id != ""):
                    continue
                if canonical_tier == 99:
                    if tier < rep_tier:
                        rep_tier = tier
                        self._rep_result = fusion_result
        fues_pep_len = {}
        for fuse_pep in fuse_peps:
            fues_pep_len[fuse_pep] = len(fuse_pep)
        sorted_pep_len = dict(sorted(fues_pep_len.items(), key=lambda item: item[1], reverse=True))
        fuse_pep_list = []
        i = 1
        for fuse_pep in sorted_pep_len:
            if fuse_pep in fused_domains:
                domains = fused_domains[fuse_pep]
            else:
                domains = predictPfam(fuse_pep, self._pfam_file)  
            types = []
            frs = fuse_peps[fuse_pep]
            for fr in frs:
                types.append(fr["type"])
            sep = ","
            type= sep.join(set(types))
            self._fuse_peps["Protein " + str(i)] = {"type": type, "length": sorted_pep_len[fuse_pep], "sequence": fuse_pep, "domains":domains, "transcripts":fuse_peps[fuse_pep]}
            i = i + 1
        
        return self._rep_result, self._fuse_peps, self._left_results, self._right_results, self._left_trans_info, self._right_trans_info 
            

    def predictPfam(self, fuse_pep):
        tmp_in = os.path.basename(tempfile.NamedTemporaryFile().name)
        out_data = ">fuse\n" + fuse_pep
        f = open(tmp_in, "w")        
        f.write(out_data)
        f.flush()
        f.close()
        script_dir = os.path.dirname(os.path.abspath(__file__))
        #os.system("hmmscan --domtblout " + tmp_out + " -E 1e-5 --domE 1e-5 " + self._pfam_file + " " + tmp_in)
        cmd = "perl " + script_dir + "/PfamScan/pfam_scan.pl -fasta " + tmp_in + " -dir " + self._pfam_file
        #print(cmd)
        stream = os.popen(cmd)
        domains = []
        for line in stream.readlines():
            line = line.strip()
            if line[0:1] != "#":
                tokens = line.split()
                if len(tokens) > 5:
                    domains.append([tokens[1], tokens[2], tokens[5], tokens[6]])
        f.close()
        os.remove(tmp_in)
        return domains    

def makeOutputHTML(in_file, out_file, temp_file, cyto_file):
    if_temp = open(temp_file, 'r')
    of_html = open(out_file,"w")
    lines = if_temp.readlines()
    sep=":"
    for line in lines:
        if line.rstrip("\n") == "//##json_data":
            in_text = open(in_file,"r")
            fusion_lines = in_text.readlines()
            json_data = {"cols":[{"title":"Details"}], "data":[]}
            details = {}
            for fusion_line in fusion_lines:
                fusion_line = fusion_line.rstrip("\n")
                fields = fusion_line.split("\t")
                if len(json_data["cols"]) == 1:
                    for field in fields:
                        if len(json_data["cols"]) < len(fields)-2:
                            title = field.replace("_", " ").title().replace("Id", "ID")
                            json_data["cols"].append({"title":title})
                else:
                    json_data["data"].append([""] + fields[:-3])
                    key = sep.join(map(str, fields[0:7]))
                    details[key] = "{\"fusion_proteins\":" + fields[-3] + ", \"left_info\":" + fields[-2] + ", \"right_info\":" + fields[-1] + "}"
            in_cytoband = open(cyto_file,"r")
            cyto_lines = in_cytoband.readlines()
            cyto_data = []
            for cyto_line in cyto_lines:
                cyto_line = cyto_line.rstrip("\n")
                cyto_fields = cyto_line.split("\t")
                if cyto_fields[0] != "chromosome":
                    cyto_data.append({"chromosome":cyto_fields[0], "bp_start":cyto_fields[1], "bp_stop":cyto_fields[2], "band":cyto_fields[3], "stain":cyto_fields[4]})
            of_html.writelines("\tvar json_data=" + json.dumps(json_data) + "\n")
            of_html.writelines("\tvar detail_data=" + json.dumps(details) + "\n")
            of_html.writelines("\tvar cyto_data=" + json.dumps(cyto_data) + "\n")
            in_text.close
        else:
            of_html.writelines(line)
    of_html.close
    if_temp.close