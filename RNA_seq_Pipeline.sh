#!/bin/bash

# Define directory paths
RAW_DIR=./rawdata/fastq
CLEAN_DIR=./cleandata
SAM_DIR=./samfile
BAM_DIR=./bamfile
SORT_BAM_DIR=./sortbam
HTML_DIR=./html
GENOME_INDEX=./Anas_platyrhynchos.ASM874695v1_hisat2_index
GTF_FILE=./Anas_platyrhynchos.ASM874695v1.109.gtf

###================== 1. Quality Control ===================
mkdir -p ${CLEAN_DIR} ${HTML_DIR}
for i in $(cat ./name.txt); do
    fastp -i ${RAW_DIR}/${i}_R1.fq.gz \
          -I ${RAW_DIR}/${i}_R2.fq.gz \
          -w 8 \
          -o ${CLEAN_DIR}/${i}_R1.clean.fq.gz \
          -O ${CLEAN_DIR}/${i}_R2.clean.fq.gz \
          -h ${HTML_DIR}/${i}.html
done

###================== 2. Alignment ===================
# Note: Genome indexing command (commented out for reference):
# hisat2-build -p 16 ./Anas_platyrhynchos.ASM874695v1.dna.toplevel.fa Anas_platyrhynchos.ASM874695v1_index

mkdir -p ${SAM_DIR}
for i in $(cat ./name.txt); do
    hisat2 -p 8 \
           -x ${GENOME_INDEX} \
           -1 ${CLEAN_DIR}/${i}_R1.clean.fq.gz \
           -2 ${CLEAN_DIR}/${i}_R2.clean.fq.gz \
           -S ${SAM_DIR}/${i}.sam
done

###================== 3. SAM to BAM Conversion ===================
# SAM to BAM conversion
mkdir -p ${BAM_DIR} ${SORT_BAM_DIR}
for i in $(cat ./name.txt); do
    samtools view -@ 8 \
                  -S ${SAM_DIR}/${i}.sam \
                  -b \
                  -o ${BAM_DIR}/${i}.bam
    samtools sort -@ 8 \
                  -o ${SORT_BAM_DIR}/${i}.sort.bam \
                  ${BAM_DIR}/${i}.bam
done
COMMENT

###================== 4. Quantification ===================
for i in $(cat ./name.txt); do
    htseq-count -f bam \
                -s no \
                ${SORT_BAM_DIR}/${i}.sort.bam \
                ${GTF_FILE} > ${i}_count.txt
done

