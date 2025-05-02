#!/bin/bash
#SBATCH -J compartment
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -c 8
#SBATCH -A csd788
#SBATCH --mem=8G
#SBATCH -t 4-00:00:00                 

# module purge
# module load cpu slurm gcc 
# export OMP_NUM_THREADS=4 
# ./pi_openmp

###############################
current="/tscc/lustre/ddn/scratch/y2xie/87.FNIH_DHC_IGM_240925"
out_dir=${current}/impute/100K/comp
raw_cell_table=${current}/all_cell_rmbl.txt
cell_table=${current}/impute/100K/impute_100k.txt

if [ -z ${genome+x} ]; then echo "genome is not specified. Default to mm10"; genome="hg38"; fi
if [[ "${genome}" == "mm10" ]]
then
    fasta_path=/tscc/projects/ps-renlab/share/bwa_indices/mm10.fa
    chrom_size=/tscc/projects/ps-renlab/y2xie/projects/genome_ref/mm10.main.chrom.sizes
    cpg=/tscc/projects/ps-renlab2/y2xie/packages/scHiCluster/cpg/mm10_cpg_ratio_100k.hdf
elif [[ "${genome}" == "hg38" ]]
then
    fasta_path=/tscc/projects/ps-renlab/share/bwa_indices/hg38.fa
    chrom_size=/tscc/projects/ps-renlab2/y2xie/projects/genome_ref/hg38.main.chrom.sizes
    cpg=/tscc/projects/ps-renlab2/y2xie/packages/scHiCluster/cpg/hg38_cpg_ratio_100k.hdf
fi
###############################

mkdir -p ${out_dir}

## Using raw matrices
/tscc/nfs/home/y2xie/anaconda3/envs/schicluster/bin/hicluster compartment --cell_table_path ${raw_cell_table} --output_prefix ${out_dir}/raw --cpg_profile_path ${cpg} --cpu 32 --resolution 10000 --chr1 1 --pos1 2 --chr2 3 --pos2 4 --chrom_size_path ${chrom_size} --mode tsv

## Using imputed matrices
/tscc/nfs/home/y2xie/anaconda3/envs/schicluster/bin/hicluster compartment --cell_table_path ${cell_table} --output_prefix ${out_dir}/impute --cpg_profile_path ${cpg} --cpu 32
