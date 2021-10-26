#!/bin/sh
#SBATCH --partition=norm,ccr
#SBATCH --cpus-per-task=64
#SBATCH --mem=80G
#SBATCH --time=24:00:00

module load python/3.7
module load hmmer/3.3.2

path=$1
patient_id=$2
case_id=$3
out_dir=$4
threads=$5

if [ -z $threads ];then
	threads="8"
fi

mkdir -p $out_dir/$patient_id/$case_id/$patient_id/db
out_file="$out_dir/$patient_id/$case_id/$patient_id/db/${patient_id}.fusion.txt"
rm -f $out_file
for fn in $path/$patient_id/$case_id/*/fusion/*.actionable.fusion.txt;do
	if [ "$fn" == "$path/$patient_id/$case_id/*/fusion/*.actionable.fusion.txt" ];then
		if [ -f "$path/$patient_id/$case_id/Actionable/${patient_id}.fusion.actionable.txt" ];then
			python fusionTools.py -i "$path/$patient_id/$case_id/Actionable/${patient_id}.fusion.actionable.txt" -o $out_file -t $threads
		else
			echo "no fusion file found"
		fi
	else
		dn=$(dirname "$fn")
		dn=$(dirname "$dn")
		sample_id=$(basename "$dn")
		exp_file=$path/$patient_id/$case_id/$sample_id/RSEM_ENS/${sample_id}.rsem_ENS.isoforms.results
		
		mkdir -p $out_dir/$patient_id/$case_id/$sample_id/fusion
		out_smp_file="$out_dir/$patient_id/$case_id/$sample_id/fusion/${sample_id}.annotated.txt"
		if [ -f $exp_file ];then
			python fusionTools.py -i $fn -m $exp_file -o $out_smp_file -t $threads
		else
			python fusionTools.py -i $fn -o $out_smp_file -t $threads
		fi
		
		if [ -f $out_file ];then
			grep -v ^left_gene $out_smp_file >> $out_file
		else
			cp $out_smp_file $out_file
		fi
	fi
done
