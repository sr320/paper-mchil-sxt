#!/bin/sh

KMER=2

LD_LIBRARY_PATH=/share/bioinformatics/gcc/lib64:/share/bioinformatics/gcc/lib

TRIN=/data/ggoetz/trinityrnaseq-2.0.6
SAMTOOLS=/share/bioinformatics/samtools/bin
BOWTIE=/share/bioinformatics/bowtie

PATH=${SAMTOOLS}:${BOWTIE}:${PATH}

${TRIN}/Trinity \
	--seqType fq \
	--max_memory 400G \
	--left mchil-sxt_1.fastq \
	--right mchil-sxt_2.fastq \
	--CPU 8 \
	--bflyCalculateCPU \
	--bflyHeapSpaceMax 16G \
	--bflyHeapSpaceInit 8G \
	--normalize_reads \
	--min_kmer_cov ${KMER} \
    --trimmomatic \
	--output trinity_All.kmer${KMER}.w_trimming \
	> trinity_output.All.kmer${KMER}.w_trimming.log
