#!/bin/bash
#$ -q shortP
#$ -j y
#$ -cwd
#$ -l h_vmem=4G
#$ -l virtual_free=4G
#$ -pe smp 12
#$ -N assem
. $HOME/.bashrc
export TEMPDIR=/genefs/MikheyevU/temp
export TMPDIR=/genefs/MikheyevU/temp
export TEMP=/genefs/MikheyevU/temp
export TMP=/genefs/MikheyevU/temp

module load trinity/r20140717 bowtie/1.1.0
Trinity --seqType fq --JM $((NSLOTS*4-5))G --single ../data/reference/raw/trimmed.fq.gz \
    --CPU $NSLOTS --output ../data/assembly/trinity
