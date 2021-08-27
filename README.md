# fusionTools
Python scripts processing fusion breakpoints and visualize them with D3.js

The objectives of fusionTools are:
1. Determine the fusion cDNA and protein sequences
2. Determine the fusion type (in-frame, out-of-frame or right gene intact)
3. Tier the importance of the fusion events
4. Visualize the results in html format

## Required packages

### Python packages

1. gtfparse
2. Bio.Seq
3. pyfaidx

### Perl packages
1. PfamScan (https://www.ebi.ac.uk/Tools/pfa/pfamscan/)

### hmmer 3.x

http://hmmer.org/download.html

## Run fusionTools

```
python fusionTools.py 
usage: fusionTools.py [-h] 
                      [--gtf GTF file] 
                      [--fasta Genome FASTA file]
                      [--canonical_trans_file Canonical transcript list]
                      [--fusion_cancer_gene_list Fusion cancer gene pair list]
                      [--cancer_gene_list Cancer gene list]
                      [--domain_file Pfam domain file]
                      [--input] --input Fusion file
                      [--output] output file 
                      [--threads Number of threads]
```

    gtf: GTF file
    fasta: Genome FASTA file
    canonical_trans_file: two column canonical transcript file (default: Refseq Mann list)
    fusion_cancer_gene_list: two column fusion gene pair list (default: Sanger Mitelman list)
    cancer_gene_list: one column cancer gene symbol list
    pfam_file: Pfam domain file
    input: input fusion list
    output: output files
    threads: number of threads
    
## Example
```
module load python
module load hmmer
export PERL5LIB=/your_pfam_scan_path:$PERL5LIB
python fusionTools -g hg19.refseq.gtf -f genome.fa -i fusion_list.txt -o processed_fusion_list.txt -t 16
```

Input example
|LeftGene|RightGene|Chr_Left|Position|Chr_Right|Position|Sample|Tool|SpanReadCount|
|------- |-------- |------- |------- |---------|--------|------|----|-------------|
|PAX7|FOXO1|chr1|19029790|chr13|41134997|RMSXXX|FusionCatcher|17|
|PAX7|FOXO1|chr1|19029790|chr13|41134997|RMSXXX|STAR-fusion|20|
|PAX7|FOXO1|chr1|19029790|chr13|41134997|RMSXXX|tophatFusion|24|
|AMD1|FARS2|chr6|111196418|chr6|5545413|RMSXXX|STAR-fusion|2|

Output example:

|left_gene|right_gene|left_chr|right_chr|left_position|right_position|sample_id|tools|type|tier|left_region|right_region|left_trans|right_trans|left_fusion_cancer_gene|right_fusion_cancer_gene|left_cancer_gene|right_cancer_gene|fusion_proteins|left_trans_info|right_trans_info|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|PAX7|FOXO1|chr1|19029790|chr13|41134997|RMS2074_D1C5FACXX|[{"FusionCatcher": 17}, {"STAR-fusion": 20}, {"tophatFusion": 24}]|in-frame|1.1|CDS|CDS|NM_001135254|NM_002015|Y|Y|Y|Y|{"MAALPGT...VSG*": {"domains": ...}}|...|...|
|AMD1|FARS2|chr6|111196418|chr6|5545413|RMS2074_D1C5FACXX|[{"STAR-fusion": 2}]|out-of-frame|4.3|CDS|CDS|NM_001634|NM_006567|N|N|N|N|{"MEAAHFF...}|...|...|

## Visaulization (not implemented yet)

![alt tag](fusionTools.png)
