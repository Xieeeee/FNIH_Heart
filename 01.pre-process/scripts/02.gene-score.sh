#!/bin/bash
#SBATCH -J genescore
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

current="/tscc/lustre/ddn/scratch/y2xie/87.FNIH_DHC_IGM_240925"
out_dir=${current}/impute/10K/genescore

gene_meta=/tscc/projects/ps-renlab/y2xie/projects/genome_ref/Paired-Tag/hg38/hg38.gcode.10X.schicluster
chrom_size=/tscc/projects/ps-renlab2/y2xie/projects/genome_ref/hg38.main.chrom.sizes

mkdir ${out_dir}
/tscc/nfs/home/y2xie/anaconda3/envs/schicluster/bin/hicluster gene-score --cell_table ${current}/impute/10K/impute_10k.txt --gene_meta ${gene_meta} --resolution 10000 --output_hdf_path ${current}/impute/10K/genescore/geneimputescore.hdf --chrom_size_path ${chrom_size} --cpu 32 --mode impute
