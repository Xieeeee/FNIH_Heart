#!/bin/bash
#SBATCH -J cellranger_arc_RNA_${start}_${end}
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 6-00:00:00                 

module purge
module load cpu slurm gcc 
# export OMP_NUM_THREADS=4 
# ./pi_openmp

### need to submit: start, end and meta
###############################
current="/tscc/projects/ps-renlab2/y2xie/projects/77.LC/81.FNIH_DPT_IGM_240827/"
# bcl2="240801_VH00454_284_AAFWJLJM5"
fastq_dir="${current}/01.rawdata/"
meta_dir="${current}/scripts/"
map_dir="${current}/03.mapping/"
script_dir="${current}/scripts/"


cellranger_mm10="/tscc/projects/ps-renlab/y2xie/projects/genome_ref/mm10"
cellranger_hg38="/tscc/projects/ps-renlab/y2xie/projects/genome_ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
cellranger_mix="/tscc/projects/ps-renlab/y2xie/projects/genome_ref/GRCh38_and_mm10"
################################

cd ${current}
mkdir -p ${map_dir}

### create the sub-meta for running
cd ${fastq_dir}
awk -v start=$start -v end=$end 'NR>=start && NR<end {print}' ${meta_dir}/${meta} > ${script_dir}/Multiome_RNA_ref_${start}_${end}.txt

while read s genome c
do
	if [[ "${genome}" == "mm10" ]]
	then 
		cellranger_ref=${cellranger_mm10}
	elif [[ "${genome}" == "hg38" ]]
	then
		cellranger_ref=${cellranger_hg38}
	elif [[ "${genome}" == "mix" ]]
        then
                cellranger_ref=${cellranger_mix}
	fi			
	/tscc/projects/ps-renlab/y2xie/packages/cellranger-6.1.2/cellranger count --id=${s} --transcriptome=${cellranger_ref} --fastqs=${fastq_dir} --sample=${s} --include-introns --chemistry=ARC-v1 #--project=${bcl2##*_} 
done <${script_dir}/Multiome_RNA_ref_${start}_${end}.txt
