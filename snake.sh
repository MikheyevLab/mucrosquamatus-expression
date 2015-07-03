#!/bin/bash
#SBATCH --job-name=align
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
##SBATCH --mail-user=%u@oist.jp
##SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=%j.out
##SBATCH --error=job_%j.err

. $HOME/.bashrc 
. ~/sasha_env/bin/activate

snakemake -j 999 -p --cluster-config cluster.json --cluster "sbatch  -p {cluster.partition} -n {cluster.n}" 