#!/bin/bash
#SBATCH -J filter
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

current="/tscc/lustre/ddn/scratch/y2xie/92.FNIH_DHC_IGM_241121"
raw_matrix=${current}/raw/single_cell
out_dir=${current}/raw/rmbl
chrom_size=/tscc/projects/ps-renlab2/y2xie/projects/genome_ref/hg38.main.chrom.sizes
bl=/tscc/projects/ps-renlab/y2xie/projects/genome_ref/hg38-blacklist.v2.bed

head -n "${end}" ${current}/all_cell_raw.txt | tail -n +"${start}" >  ${current}/tmp_${start}_${end}.txt
mkdir -p ${out_dir}
/tscc/nfs/home/y2xie/anaconda3/envs/schicluster/bin/hicluster filter-contact --output_dir ${out_dir} --blacklist_1d_path ${bl} --chr1 1 --pos1 2 --chr2 3 --pos2 4 --contact_table ${current}/tmp_${start}_${end}.txt --chrom_size_path ${chrom_size} --not_remove_duplicates
