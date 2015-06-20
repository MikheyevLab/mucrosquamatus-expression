#!/bin/bash
#$ -q genomics
#$ -j y
#$ -cwd
##$ -pe smp 16
#$ -l h_vmem=400G
#$ -l virtual_free=400G
#$ -N nb

. $HOME/.bashrc
export TEMPDIR=/genefs/MikheyevU/temp
export TMPDIR=/genefs/MikheyevU/temp
export TEMP=/genefs/MikheyevU/temp
export TMP=/genefs/MikheyevU/temp

base=Pm

#newAssembly -cdna -force $base
cd $base
#addRun -lib illumina /genefs/MikheyevU/sasha/projects/mucrosquamatus-expression/data/reference/raw/trimmed.fq
runProject -m -cpu 16 -mi 98 -siom 390 -a 200 -urt -novs -het

#resuming on genomics after Newbler ran out of memory

#cd assembly
#runProject -large -m -cpu 1 -mi 95 -siom 240 

#-novs no vector screening
#-l minimum length in large contigs
#-mi minimum % identity
#-m everything in memory
#-ud no duplicates
#-a minimum contig length
#-nohet no heterozygosity
#-novs no vector screening
#-urt use read tips
#-het/-nohet  set heterozygosity
#-minlen minimum read length up to 45
