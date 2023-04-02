#!/bin/bash

# chimpanzee_align
# 1.0.0

fq=$1
experimenttype=$2
referencegenome=$3
referencete=$4
minimap2=$5

# align fastq reads to reference genome and te genome using minimap2
# use different parameters depending on experiment type (DCS or DRS)
# sort alignments using samtools
if [[ $experimenttype == 'DCS' ]]
then
  prefix=$(basename $fq | sed 's/.fastq//')
  $minimap2 -ax splice --secondary=no --sam-hit-only -Y $referencegenome $fq > ./alignment/genome/aligned_genome_${prefix}.sam
  $minimap2 -ax splice --secondary=no --sam-hit-only -Y $referencete $fq > ./alignment/te/aligned_te_${prefix}.sam
  module load samtools
  samtools sort -T ${prefix}_gtemp -o ./alignment/genome/sorted_aligned_genome_${prefix}.sam ./alignment/genome/aligned_genome_${prefix}.sam
  wait
  samtools sort -T ${prefix}_ttemp -o ./alignment/te/sorted_aligned_te_${prefix}.sam ./alignment/te/aligned_te_${prefix}.sam
  wait
elif [[ $experimenttype == 'DRS' ]]
then
  prefix=$(basename $fq | sed 's/.fastq//')
  $minimap2 -ax splice -uf -k14 --secondary=no --sam-hit-only -Y $referencegenome $fq > ./alignment/genome/aligned_genome_${prefix}.sam
  $minimap2 -ax splice -uf -k14 --secondary=no --sam-hit-only -Y $referencete $fq > ./alignment/te/aligned_te_${prefix}.sam
  module load samtools
  samtools sort -T ${prefix}_gtemp -o ./alignment/genome/sorted_aligned_genome_${prefix}.sam ./alignment/genome/aligned_genome_${prefix}.sam
  wait
  samtools sort -T ${prefix}_ttemp -o ./alignment/te/sorted_aligned_te_${prefix}.sam ./alignment/te/aligned_te_${prefix}.sam
  wait
fi
