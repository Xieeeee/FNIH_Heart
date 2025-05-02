#!/bin/bash
#SBATCH -J split_hic_${s}
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -N 1
#SBATCH -A csd788
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 6-00:00:00

module purge
module load cpu slurm gcc
# export OMP_NUM_THREADS=4
# ./pi_openmp

#################################
inpath="/tscc/projects/ps-renlab2/y2xie/projects/77.LC/87.FNIH_DHC_IGM_240925/"
subdir="./"
#################################

if [ -z ${mode+x} ]; then echo "10X kit used is not specified. Default to arc"; mode="arc"; fi
bash /tscc/projects/ps-renlab2/y2xie/scripts/Paired-HiC/03.preproc_paired_hic_tscc2.sh -i ${s} -p ${inpath} -d ${subdir} -m ${mode}