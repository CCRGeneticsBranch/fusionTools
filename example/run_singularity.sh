singularity exec --bind /data/khanlab/projects/processed_DATA,/data/Clinomics/Ref/khanlab/,$PWD ../fusion_tools_v1.sif processFusionCase.sh \
                       -d /data/khanlab/projects/processed_DATA \
                       -p RH4 \
                       -c RNA_RMS_HiBiT \
                       -o $PWD/out_dir2 \
                       -f /data/Clinomics/Ref/khanlab/PfamDB \
                       -g /data/Clinomics/Ref/khanlab/ucsc.hg19.fasta

