#!/bin/bash
#SBATCH -J split
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -A csd788
#SBATCH --mem=8G
#SBATCH -t 4-00:00:00                 

# module purge
# module load cpu slurm gcc 
# export OMP_NUM_THREADS=4 
# ./pi_openmp

current="/tscc/projects/ps-renlab2/y2xie/projects/77.LC/87.FNIH_DHC_IGM_240925"
bash /tscc/projects/ps-renlab2/y2xie/scripts/Paired-HiC/phc.batch_splitPairs_single.sh -i ${current}/03.mapping -o /tscc/lustre/ddn/scratch/y2xie/87.FNIH_DHC_IGM_240925/single_cell/ -m ${current}/05.R/FNIH_Liver_DHC_demultiplex_SNG.cluster_info_new.txt -s ${start} -e ${end}
