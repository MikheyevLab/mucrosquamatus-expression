#!/bin/bash
#$ -q shortP
#$ -j y
#$ -cwd
#$ -l h_vmem=4G
#$ -l virtual_free=4G
#$ -pe smp 12
#$ -N tt
. $HOME/.bashrc
export TEMPDIR=/genefs/MikheyevU/temp
export TMPDIR=/genefs/MikheyevU/temp
export TEMP=/genefs/MikheyevU/temp
export TMP=/genefs/MikheyevU/temp

module load trinity/r20140717 bowtie/1.1.0
TRINITY_HOME=/apps/gnu/trinity/r20140717/bin

$TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts ../data/assembly/trinity/Trinity.fasta --seqType fq --single ../data/reference/raw/trimmed.fq --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference --thread_count $NSLOTS
