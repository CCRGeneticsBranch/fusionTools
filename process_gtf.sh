#!/usr/bin/bash

#./process_gtf.sh data/gencode.v36lift37.annotation.gtf.gz data/gencode.v36lift37.pc_translations.fa.gz
gtf=$1
pfa=$2
prefix="${gtf%.gtf.gz}"
prefixp="${pfa%.fa.gz}"
if [ ! -f ${prefix}.sorted.gtf.gz ];then
	zcat $gtf | grep -v ^"#" | sort -k1,1 -k4,4n > ${prefix}.sorted.gtf
	#perl gen_canonical_gtf.pl -s ${prefix}.sorted.gtf -c ${prefix}.ensembl_canonical_mane.txt -o ${prefix}.sorted.genename_changed.gtf -g ${prefix}.sorted.genename_changed.canonical.gtf
	module load samtools
	module load python
	bgzip < ${prefix}.sorted.gtf > ${prefix}.sorted.gtf.gz
	#bgzip < ${prefix}.sorted.genename_changed.gtf > ${prefix}.sorted.genename_changed.gtf.gz
	#bgzip < ${prefix}.sorted.genename_changed.canonical.gtf > ${prefix}.sorted.genename_changed.canonical.gtf.gz
	tabix -p gff ${prefix}.sorted.gtf.gz
	#tabix -p gff ${prefix}.sorted.genename_changed.gtf.gz
	#tabix -p gff ${prefix}.sorted.genename_changed.canonical.gtf.gz
fi
module load python/3.7
mkdir -p tmp
rm -f tmp/*
python convertGeneBED.py --gtf ${prefix}.sorted.gtf.gz --bed ${prefix}.sorted.genes.bed
python makeDomainFile.py --in_file $pfa --out_file $prefixp.tsv.tmp
python saveTranscriptTable.py --gtf ${prefix}.sorted.gtf.gz --out ${prefix}.transcripts.tsv
split -l 3000 $prefixp.tsv.tmp tmp/gencode_
mkdir -p tmp/out
for fn in tmp/gencode*;do
	dn=$(dirname "$fn")
	bn=$(basename "$fn")
	sbatch makeDomainFileChunk.sh $fn tmp/out/${bn}.txt
done
#echo -e "trans\tprotein_length\tdomains" > data/gencode.v36lift37.domains.tsv
#cat tmp/out/*.txt >> data/gencode.v36lift37.domains.tsv



