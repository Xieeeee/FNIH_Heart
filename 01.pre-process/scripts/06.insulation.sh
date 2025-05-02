#!/bin/bash
#SBATCH -J insulation
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -c 8
#SBATCH -A csd788
#SBATCH --mem=8G
#SBATCH -t 4-00:00:00                 

source /tscc/nfs/home/y2xie/anaconda3/etc/profile.d/conda.sh
conda activate schicluster

###############################
current="/tscc/lustre/ddn/scratch/y2xie/87.FNIH_DHC_IGM_240925"
out_dir=${current}/impute/25K/domain
cell_table=${current}/impute/25K/impute_25k.txt
###############################

mkdir -p ${out_dir}
hicluster domain --cell_table_path ${cell_table} --output_prefix ${out_dir}/domain --resolution 25000 --window_size 10 --cpu 32
