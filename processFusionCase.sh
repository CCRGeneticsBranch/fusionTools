#!/bin/bash
#SBATCH --partition=norm,ccr
#SBATCH --cpus-per-task=64
#SBATCH --mem=80G
#SBATCH --time=24:00:00

while [[ "$#" -gt 0 ]]; do case $1 in
  -d) path="$2"; shift;;
  -p) patient_id="$2"; shift;;
  -c) case_id="$2"; shift;;
  -f) pfam_path="$2"; shift;;
  -g) genome_fasta="$2"; shift;;
  -t) threads="$2"; shift;;
  -o) out_dir="$2"; shift;;
  -v) gencode_version="$2"; shift;;
  *) echo "Unknown parameter passed: $1"; exit 1;;
esac; shift; done

if [ -z $path ] || [ -z $patient_id ] || [ -z $case_id ] || [ -z $pfam_path ] || [ -z $genome_fasta ];then
	echo "usage: processFusionCase.h "
	echo
	echo "required:"
	echo "-d: processed data path"
	echo "-p: patient ID"
	echo "-c: case ID"
	echo "-f: Pfam DB folder"
	echo "-g: Genome fasta"
	echo 
	echo "optional:"
	echo "-t: number of threads. (default: SLURM_CPUS_PER_TASK variable)"
	echo "-o: output folder (default: same as input folder)"
	echo "-v: Gencode version (default: 37)"
	echo
	echo "Example:"
	echo
	echo "./processFusionCase.sh -d /data/Compass/Analysis/ProcessedResults_NexSeq/ExomeRNA_Results \\"
	echo "                       -p CP02796 \\"
	echo "                       -c RT-0391 \\"
	echo "                       -f /data/Clinomics/Ref/khanlab/PfamDB \\"
	echo "                       -g /data/Clinomics/Ref/khanlab/ucsc.hg19.fasta \\"
	echo "                       -v 36"
	echo 
	exit 0
fi

script_file=`realpath $0`
script_home=`dirname $script_file`

if [ -z $gencode_version ];then
	gencode_version="37"
fi

gtf=${script_home}/data/gencode.v${gencode_version}lift37.annotation.sorted.gtf.gz
can_file=${script_home}/data/gencode.v${gencode_version}lift37.canonical.txt
domain_file=${script_home}/data/gencode.v${gencode_version}lift37.domains.tsv
if [ ! -f $gtf ];then
	echo "GTF $gtf does not exist"
	exit 1
fi
if [ ! -f $can_file ];then
	echo "Canonical file $gtf does not exist"
	exit 1
fi
if [ ! -f $domain_file ];then
	echo "Domain file $domain_file does not exist"
	exit 1
fi

if [ -z $threads ];then
	threads=$SLURM_CPUS_PER_TASK
	if [ -z $threads ];then
		threads=4
	fi
fi

if [ -z $out_dir ];then
	out_dir=$path
fi

mkdir -p $out_dir/$patient_id/$case_id/$patient_id/db
out_file="$out_dir/$patient_id/$case_id/$patient_id/db/${patient_id}.fusion"
rm -f $out_file
for fn in $path/$patient_id/$case_id/*/fusion/*.actionable.fusion.txt;do
	if [ "$fn" == "$path/$patient_id/$case_id/*/fusion/*.actionable.fusion.txt" ];then
		if [ -f "$path/$patient_id/$case_id/Actionable/${patient_id}.fusion.actionable.txt" ];then
			python $script_home/fusionTools.py -i "$path/$patient_id/$case_id/Actionable/${patient_id}.fusion.actionable.txt" -o $out_file -t $threads -p $pfam_path -f $genome_fasta -g $gtf -n $can_file -d $domain_file
		else
			echo "no fusion file found"
		fi
	else
		dn=$(dirname "$fn")
		dn=$(dirname "$dn")
		sample_id=$(basename "$dn")
		exp_file=$path/$patient_id/$case_id/$sample_id/RSEM_ENS/${sample_id}.rsem_ENS.isoforms.results
		if [ ! -f $exp_file ];then
			exp_file=$path/$patient_id/$case_id/$sample_id/RSEM/${sample_id}.rsem.isoforms.results
		fi
		mkdir -p $out_dir/$patient_id/$case_id/$sample_id/fusion
		out_smp_file="$out_dir/$patient_id/$case_id/$sample_id/fusion/${sample_id}.annotated"
		if [ -f $exp_file ];then
			python $script_home/fusionTools.py -i $fn -m $exp_file -o $out_smp_file -t $threads -p $pfam_path -f $genome_fasta  -g $gtf -n $can_file -d $domain_file
		else 
			python $script_home/fusionTools.py -i $fn -o $out_smp_file -t $threads -p $pfam_path -f $genome_fasta -g $gtf -n $can_file -d $domain_file
		fi
		
		if [ -f ${out_file}.txt ];then
			grep -v ^left_gene ${out_smp_file}.txt >> ${out_file}.txt
		else
			cp ${out_smp_file}.txt ${out_file}.txt
		fi		
	fi
done
python $script_home/makeOutputHTML.py -i ${out_file}.txt -o ${out_file}.html -t $script_home/data/template.html -c $script_home/data/hg19_cytoBand.txt