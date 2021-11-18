perl gen_canonical_gtf.pl -s gencode.v38lift37.annotation.sorted.gtf -c data/Ensembl_canonical_mane.txt -o gencode.v38lift37.annotation.sorted.genename_changed.gtf -g gencode.v38lift37.annotation.sorted.genename_changed.canonical.gtf

bgzip < gencode.v38lift37.annotation.sorted.genename_changed.gtf > gencode.v38lift37.annotation.sorted.genename_changed.gtf.gz
bgzip < gencode.v38lift37.annotation.sorted.genename_changed.canonical.gtf > gencode.v38lift37.annotation.sorted.genename_changed.canonical.gtf.gz
module load samtools
tabix -p gff gencode.v38lift37.annotation.sorted.genename_changed.gtf.gz
tabix -p gff gencode.v38lift37.annotation.sorted.genename_changed.canonical.gtf.gz