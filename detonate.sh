#!/bin/bash
#$ -q longP
#$ -j y
#$ -cwd
#$ -N td
#$ -pe smp 12
#$ -l h_vmem=2G
#$ -l virtual_free=2G

. $HOME/.bashrc

module load bowtie/1.1.0
assembly=../data/assembly/trinity/Trinity.fasta ; out=trinity
#assembly=../data/assembly/newbler/assembly/454AllContigs.fna; out=newbler
/apps/MikheyevU/sasha/detonate-1.9/rsem-eval/rsem-eval-calculate-score  ../data/reference/raw/trimmed.fq $assembly $out 143 -p $NSLOTS
