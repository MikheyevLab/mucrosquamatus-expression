#!/bin/bash
#$ -q shortP
#$ -j y
#$ -cwd
#$ -l h_vmem=1G
#$ -l virtual_free=1G
#$ -pe smp 16 
#$ -N tr
. $HOME/.bashrc
export TEMPDIR=/genefs/MikheyevU/temp
export TMPDIR=/genefs/MikheyevU/temp
export TEMP=/genefs/MikheyevU/temp
export TMP=/genefs/MikheyevU/temp

bbmap=/imports/hpc_software/MikheyevU/sasha/bbmap 

$bbmap/bbduk2.sh in=../data/reference/raw/snake_S1_L001_R1_001.fastq.gz out=../data/reference/raw/trimmed.fq.gz ref=$bbmap/resources/nextera.fa.gz \
lliteral=GTTTTCCCAGTCACGACAATTGCAGTGGTATCAACGCAGAGCGGCCGCTTTT,AGTGGAGAGCTAACAATTTCACACAGGAAAGCAGTGGTATCAACGCAGAGTACATGG,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA \
threads=$NSLOTS \
mink=12 \
hdist=1 \
qtrim=r trimq=10 \
-Xmx"$(($NSLOTS-5))"G

