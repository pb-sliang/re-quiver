#!/bin/bash
#$ -N aln-quiver
#$ -cwd
#$ -j y
#$ -q bigmem
#$ -pe smp 24
#$ -l h_rt=24:0:0,h_vmem=50G,mem_free=15G
#$ -o log
#$ -e log


source /mnt/software/Modules/current/init/bash

module load smrtanalysis/mainline
module load samtools

set -vex

pbalign --tmpDir=scratch/ --nproc=16 --minAccuracy=0.75 --minLength=50 --minAnchorSize=12 --maxDivergence=30 --concordant --algorithm=blasr          --algorithmOptions=--useQuality --maxHits=1 --hitPolicy=randombest --seed 0 ${ctg}.bam ${ctg}_ref.fa scratch/aln-${ctg}.bam

samtools faidx ${ctg}_ref.fa

variantCaller --algorithm=arrow -x 5 -X 120 -q 20 -j 16 -r ${ctg}_ref.fa scratch/aln-${ctg}.bam  -o cns-${ctg}.fasta.gz -o cns-${ctg}.fastq.gz
