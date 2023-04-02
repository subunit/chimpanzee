#!/bin/bash

# chimpanzee
# 1.0.0

# sample information
samplename=
experimenttype=
species=

# input files
inputdir=
referencegenome=
referencete=
geneanno=
refseqid=

# software
minimap2=
runpychopper=

# miscellaneous
netid=
minimumhit=

# load dependencies
echo -e "$(date)\t\tloading modules"
module purge
module load samtools
module load bedtools2
module load R
wait
echo -e "$(date)\t\tmodules loaded\n"

# set environment
cd $inputdir
mkdir chimpanzee
mkdir fastq_input
mkdir misc
mv *.fastq ./fastq_input
fastq_input=./fastq_input
wait

# run pychopper
if [[ $experimenttype == 'DCS' ]] && [[ $runpychopper == 'Y' ]]
then
  echo -e "$(date)\t\tpreprocessing and orienting reads\n"
  mkdir -p pychopper/{classified,unclassified}
  mkdir -p misc/pychopper
  uniquetime=`date +"%m%d%y%H%M"`
  for i in ./fastq_input/*.fastq
  do
    prefix=$(basename $i | sed 's/.fastq//')
    sbatch --mem=4G --cpus-per-task=1 --ntasks=1 --job-name="${uniquetime}_pychopper" -Q --output="preprocess_%A.out" --wrap="cdna_classifier.py -r ./misc/pychopper/report_${prefix}.pdf -u ./pychopper/unclassified/unclassified_${prefix}.fastq -w ./pychopper/classified/rescued_${prefix}.fastq $i ./pychopper/classified/classified_${prefix}.fastq"
  done
  sleep 60
  echo ""
  while [[ $(squeue -h -o %j -u $netid | grep -c "${uniquetime}_pychopper") != 0 ]]
  do
    echo -e "$(date)\t\tpychopper in progress ... $(squeue -h -o %j -u $netid | grep -c "${uniquetime}_pychopper") job(s) remaining"
    sleep 60
  done
  wait
  fastq_input=./pychopper/classified
  rm preprocess_*
  echo -e "$(date)\t\treads preprocessed and oriented"
fi

# align fastq reads to genome and te with minimap2
echo -e "$(date)\t\taligning reads\n"
mkdir -p alignment/{genome,te}
uniquetime=`date +"%m%d%y%H%M"`
for i in ${fastq_input}/*.fastq
do
  sbatch --mem=32G --cpus-per-task=4 --ntasks=1 --job-name="${uniquetime}_map" -Q --output="map_%A.out" chimpanzee_align.sh $i $experimenttype $referencegenome $referencete $minimap2
done
sleep 60
echo ""
while [[ $(squeue -h -o %j -u $netid | grep -c "${uniquetime}_map") != 0 ]]
do
  echo -e "$(date)\t\talignment in progress ... $(squeue -h -o %j -u $netid | grep -c "${uniquetime}_map") job(s) remaining"
  sleep 30
done
rm map_*
echo -e "$(date)\t\treads aligned\n"

# merge sam files in batches
echo -e "$(date)\t\tmerging genomic reads\n"
uniquetime=`date +"%m%d%y%H%M"`
rm ./alignment/genome/aligned_genome_*
mkdir -p alignment/genome_merge_sup
fileno=1
dirno=1
mkdir -p alignment/genome_merge_slice_${dirno}
for i in ./alignment/genome/sorted_aligned_genome_*
do
  if [[ $fileno == 51 ]]
  then
    sbatch --mem=8G --cpus-per-task=1 --ntasks=1 --job-name="${uniquetime}_gmerge" -Q --output="gmerge_%A.out" --wrap="samtools merge ./alignment/genome_merge_sup/genome_merge_sup_${dirno}.sam ./alignment/genome_merge_slice_${dirno}/*"
    dirno=$(( dirno + 1 ))
    mkdir -p alignment/genome_merge_slice_${dirno}
    fileno=1
  fi
  cp $i ./alignment/genome_merge_slice_${dirno}/
  fileno=$(( fileno + 1 ))
done
if [[ $(ls -l ./alignment/genome_merge_slice_${dirno}/*.sam | wc -l) == 1 ]]
then
  cp ./alignment/genome_merge_slice_${dirno}/*.sam ./alignment/genome_merge_sup/
else
  samtools merge ./alignment/genome_merge_sup/genome_merge_sup_${dirno}.sam ./alignment/genome_merge_slice_${dirno}/*
fi
echo ""
while [[ $(squeue -h -o %j -u $netid | grep -c "${uniquetime}_gmerge") != 0 ]]
do
  echo -e "$(date)\t\tmerging in progress ... $(squeue -h -o %j -u $netid | grep -c "${uniquetime}_gmerge") job(s) remaining"
  sleep 30
done
rm gmerge_*
samtools merge ./chimpanzee/${samplename}_aligned_genome.sam ./alignment/genome_merge_sup/*
wait
rm -r ./alignment/genome_merge_*
echo -e "$(date)\t\tgenomic reads merged\n"

echo -e "$(date)\t\tmerging te reads\n"
uniquetime=`date +"%m%d%y%H%M"`
rm ./alignment/te/aligned_te_*
mkdir -p alignment/te_merge_sup
fileno=1
dirno=1
mkdir ./alignment/te_merge_slice_${dirno}
for i in ./alignment/te/sorted_aligned_te_*
do
  if [[ $fileno == 201 ]]
  then
    sbatch --mem=8G --cpus-per-task=1 --ntasks=1 --job-name="${uniquetime}_tmerge" -Q --output="tmerge_%A.out" --wrap="samtools merge ./alignment/te_merge_sup/te_merge_sup_${dirno}.sam ./alignment/te_merge_slice_${dirno}/*"
    dirno=$(( dirno + 1 ))
    mkdir -p alignment/te_merge_slice_${dirno}
    fileno=1
  fi
  cp $i ./alignment/te_merge_slice_${dirno}/
  fileno=$(( fileno + 1 ))
done
if [[ $(ls -l ./alignment/te_merge_slice_${dirno}/*.sam | wc -l) == 1 ]]
then
  cp ./alignment/te_merge_slice_${dirno}/*.sam ./alignment/te_merge_sup/
else
  samtools merge ./alignment/te_merge_sup/te_merge_sup_${dirno}.sam ./alignment/te_merge_slice_${dirno}/*
fi
echo ""
while [[ $(squeue -h -o %j -u $netid | grep -c "${uniquetime}_tmerge") != 0 ]]
do
  echo -e "$(date)\t\tmerging in progress ... $(squeue -h -o %j -u $netid | grep -c "${uniquetime}_tmerge") job(s) remaining"
  sleep 15
done
rm tmerge_*
samtools merge ./chimpanzee/${samplename}_aligned_te.sam ./alignment/te_merge_sup/*
wait
rm -r ./alignment/te_merge_*
echo -e "$(date)\t\tte reads merged"

# extract sam header for genomic alignment
# sort alignment by query name
echo -e "$(date)\t\tsorting genomic alignments"
samtools view -H ./chimpanzee/${samplename}_aligned_genome.sam > ./misc/${samplename}_genome_header.txt
samtools sort -T gtemp -n -o ./chimpanzee/${samplename}_aligned_genome.sorted.sam ./chimpanzee/${samplename}_aligned_genome.sam
wait
# samtools sort -T gtemp2 -o ${samplename}_aligned_genome_coor.sorted.sam ${samplename}_aligned_genome_split.sam
# wait
# grep -v '^@' ${samplename}_aligned_genome_coor.sorted.sam > ${samplename}_aligned_genome_coor_headerless.sorted.sam
echo -e "$(date)\t\tgenomic alignment sorted\n"
echo -e "$(date)\t\tsorting te alignments"
samtools sort -T ttemp -n -o ./chimpanzee/${samplename}_aligned_te.sorted.sam ./chimpanzee/${samplename}_aligned_te.sam
wait
grep -v '^@' ./chimpanzee/${samplename}_aligned_te.sorted.sam > ./misc/${samplename}_aligned_te_headerless.sorted.sam
echo -e "$(date)\t\tte alignment sorted\n"

# collapse repeat query names
# generate common query names between genomic alignment and te alignment
# remove leftover sam header
echo -e "$(date)\t\tgenerating common reads"
awk -F'\t' '{print $1}' ./chimpanzee/${samplename}_aligned_genome.sorted.sam | uniq > ./misc/${samplename}_genome_unique_hits.txt
awk -F'\t' '{print $1}' ./misc/${samplename}_aligned_te_headerless.sorted.sam | uniq > ./misc/${samplename}_te_unique_hits.txt
awk 'NR==FNR{common[$0]; next} $0 in common' ./misc/${samplename}_genome_unique_hits.txt ./misc/${samplename}_te_unique_hits.txt > ./misc/${samplename}_common_reads.txt
wait
echo -e "$(date)\t\tcommon reads generated\n"

# split genomic alignment into manageable chunks for parallel processing
echo -e "$(date)\t\tslicing sam file"
split -l 20000 -d -a 4 ./chimpanzee/${samplename}_aligned_genome.sorted.sam ./misc/samsliced
wait
echo -e "$(date)\t\tsam file sliced into $(ls -l ./misc/samsliced* | wc -l) files\n"

# match alignments with common query name list to extract common alignments
echo -e "$(date)\t\textracting common reads\n"
uniquetime=`date +"%m%d%y%H%M"`
for i in ./misc/samsliced*
do
  prefix=$(basename $i)
  sbatch --mem=8G --cpus-per-task=1 --ntasks=1 --job-name="${uniquetime}_extract" -Q --output="extract_%A.out" --wrap="grep -f ./misc/${samplename}_common_reads.txt $i > ./misc/common_${prefix}.txt"
done
sleep 60
echo ""
while [[ $(squeue -h -o %j -u $netid | grep -c "${uniquetime}_extract") != 0 ]]
do
  echo -e "$(date)\t\textraction in progress ... $(squeue -h -o %j -u $netid | grep -c "${uniquetime}_extract") job(s) remaining"
  sleep 30
done
cat ./misc/common_samsliced* > ./misc/${samplename}_common_genome_headerless.sam
wait
rm extract_*
rm ./misc/samsliced*
rm ./misc/common_samsliced*
echo -e "$(date)\t\tcommon reads extracted\n"

# final processing
echo -e "$(date)\t\tprocessing data"
# sort common alignments
cat ./misc/${samplename}_genome_header.txt ./misc/${samplename}_common_genome_headerless.sam > ./chimpanzee/${samplename}_common_genome.sam
samtools view -b ./chimpanzee/${samplename}_common_genome.sam > ./misc/${samplename}_common_genome.bam
samtools sort -T gbtemp -m 16G -o ./chimpanzee/${samplename}_common_genome.sorted.bam ./misc/${samplename}_common_genome.bam
wait
# convert sorted bam file to sorted bed file
bedtools bamtobed -bed12 -i ./chimpanzee/${samplename}_common_genome.sorted.bam > ./misc/${samplename}_common_genome.bed
# awk 'NR==FNR{arr[NR]=$10;next} {$9=arr[FNR]} 1' ${samplename}_aligned_genome_coor_headerless.sorted.sam ${samplename}_common_genome.bed > ${samplename}_common_genome_seq.bed
bedtools sort -i ./misc/${samplename}_common_genome.bed > ./chimpanzee/${samplename}_common_genome.sorted.bed
wait
# annotate sorted bed file
bedtools intersect -wa -wb -a ./chimpanzee/${samplename}_common_genome.sorted.bed -b $geneanno | cut -f 4,16 > ./misc/${samplename}_common_genome_anno.txt
# slice out duplicates in te alignments
awk '{dup[$1]++; if(dup[$1]==1){item[$1]=$0}; if(dup[$1]==2){print item[$1]}}' ./misc/${samplename}_aligned_te_headerless.sorted.sam > ./misc/${samplename}_aligned_te_dup.sam
# pair gene with te
awk 'NR==FNR{te[$1]=$3; next} {for(x in te){if($1==x){print $0 "\t" te[x]}}}' ./misc/${samplename}_aligned_te_headerless.sorted.sam ./misc/${samplename}_common_genome_anno.txt > ./misc/${samplename}_common_genome_te_anno_norep.txt
awk 'NR==FNR{te[$1]=$3; next} {for(x in te){if($1==x){print $0 "\t" te[x]}}}' ./misc/${samplename}_aligned_te_dup.sam ./misc/${samplename}_common_genome_anno.txt > ./misc/${samplename}_common_genome_te_anno_dup.txt
cat ./misc/${samplename}_common_genome_te_anno_norep.txt ./misc/${samplename}_common_genome_te_anno_dup.txt > ./misc/${samplename}_common_genome_te_anno.txt
cut -f 2,3 ./misc/${samplename}_common_genome_te_anno.txt | sort -k 1,1 -k 2,2 > ./chimpanzee/${samplename}_chimera.txt
# awk '{count[$1"\t"$2]++;if(count[$1"\t"$2]==1){first[$1"\t"$2]=$0};if(count[$1"\t"$2]==5){print first[$1"\t"$2]"\n"first[$1"\t"$2]"\n"first[$1"\t"$2]"\n"first[$1"\t"$2]};if(count[$1"\t"$2]>4){print}}' ${samplename}_common_genome_te_anno_final.sorted.txt > ${samplename}_common_genome_te_anno_final_sig.sorted.txt
# retain hits above threshold
awk -v x=$minimumhit '{rep[$0]++; if(rep[$0]==1){item[$0]=$0}; if(rep[$0]==x){for(i=1; i<=x; i++){print item[$0]}}; if(rep[$0]>x){print}}' ./chimpanzee/${samplename}_chimera.txt > ./chimpanzee/${samplename}_chimera_sig.txt
# awk 'BEGIN{x=1}{count[$0]++;if(count[$0]==1){first[$0]=$0};if(count[$0]==50){while(x<count[$0]){print first[$0];x++};x=1};if(count[$0]>49){print}}' ${samplename}_common_genome_te_anno_final.sorted.txt > ${samplename}_common_genome_te_anno_final_sig.sorted.txt
# convert refsed id to gene symbol
cut -f 1 ./chimpanzee/${samplename}_chimera_sig.txt | cut -f 1 -d '.' > ./misc/${samplename}_chimera_sig_refseq.txt
awk 'NR==FNR{conv[$1]=$2; next} {print conv[$0]}' $refseqid ./misc/${samplename}_chimera_sig_refseq.txt > ./misc/${samplename}_chimera_sig_geneid.txt
cut -f 2 ./chimpanzee/${samplename}_chimera_sig.txt | paste ./misc/${samplename}_chimera_sig_geneid.txt - > ./chimpanzee/${samplename}_chimera_annotated.txt
echo -e "$(date)\t\tdata processed\n"
