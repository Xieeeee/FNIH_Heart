#!/bin/bash
#SBATCH -J bed2bam2bw
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -N 1
#SBATCH -A csd788
#SBATCH -n 4
#SBATCH -c 4
#SBATCH --mem=4G
#SBATCH -t 12:00:00

module purge
module load cpu slurm gcc
# export OMP_NUM_THREADS=4
# ./pi_openmp

genome=hg38
if [[ "${genome}" == "mm10" ]]
then
    bl=/tscc/projects/ps-renlab2/y2xie/ps-renlab/y2xie/projects/genome_ref/mm10_CnR_blacklist.bed
    size=2652783500
elif [[ "${genome}" == "hg38" ]]
then
    bl=/tscc/projects/ps-renlab2/y2xie/ps-renlab/y2xie/projects/genome_ref/hg38_CnR_blacklist.bed
    size=2913022398
fi

### bam to tagalign: python /tscc/projects/ps-renlab2/y2xie/git/peak-call-pipeline/scripts/make_tagAligns.py -b $f -o ${fname}_
dir=/tscc/projects/ps-renlab2/y2xie/projects/77.LC/81.FNIH_DPT_IGM_240827/03.mapping/10X/ATAC_split_bam/chamber
bedToBam -i ${dir}/${s}.tagAlign.gz -g /tscc/projects/ps-renlab2/y2xie/projects/genome_ref/hg38.main.chrom.sizes > ${dir}/${s}.bam
samtools sort -@ 16 -o ${dir}/${s}_sorted.bam ${dir}/${s}.bam
samtools index ${dir}/${s}_sorted.bam
# bash /tscc//projects/ps-renlab2/y2xie/scripts/DPT/12.bamCoverage.sh -i ${dir}/${s}_sorted.bam -g hg38 -o ${dir}
bamCoverage -b ${dir}/${s}_sorted.bam -o ${dir}/${s}.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize ${size} --ignoreForNormalization chrY --extendReads -p 2 --extendReads 200 -bl ${bl} --minFragmentLength 20
rm ${dir}/${s}.bam
