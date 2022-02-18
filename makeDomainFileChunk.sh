#!/usr/bin/bash
#SBATCH --partition=norm,ccr
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
module load python/3.7
module load hmmer

python makeDomainFileChunk.py -i $1 -o $2
