#./gen_canonical_gtf.sh v39
module load samtools
version=$1
perl gen_canonical_gtf.pl -s data/gencode.${version}.annotation.sorted.gtf -c data/gencode.${version}.canonical.txt -o gencode.${version}.annotation.sorted.genename_changed.gtf -g gencode.${version}.annotation.sorted.genename_changed.canonical.gtf

bgzip < gencode.${version}.annotation.sorted.genename_changed.gtf > gencode.${version}.annotation.sorted.genename_changed.gtf.gz
bgzip < gencode.${version}.annotation.sorted.genename_changed.canonical.gtf > gencode.${version}.annotation.sorted.genename_changed.canonical.gtf.gz
module load samtools
tabix -p gff gencode.${version}.annotation.sorted.genename_changed.gtf.gz
tabix -p gff gencode.${version}.annotation.sorted.genename_changed.canonical.gtf.gz