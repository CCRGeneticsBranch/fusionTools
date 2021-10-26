#!/bin/sh

for fn in data/ensembl_domains/ens*;do
	dn=$(dirname "$fn")
	bn=$(basename "$fn")
	sbatch makeDomainFileChunk.sh $fn $dn/out_${bn}.txt	
done