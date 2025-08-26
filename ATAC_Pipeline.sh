#!/bin/bash

# Define directory paths
RAW_DIR=/storage-04/duck_library/follicle/01.RawData
CLEAN_DIR=./cleandata
BAM_DIR=./sortbam
DEDUP_DIR=./dedup_bam
PEAK_DIR=./callpeak
BW_DIR=./bw
TSS_DIR=./tss
GENOME_REF=/storage-02/lizhen/genome/Anas_platyrhynchos.ASM874695v1.dna.toplevel.fa
GENOME_BED=/storage-02/lizhen/genome/Anas_platyrhynchos.ASM874695v1.110.bed

###================== 1. Quality Control ===================
# Initial QC (commented out)
# fastqc -o ./fastqc ./rawdata/* -t 5

# Read trimming with fastp
mkdir -p ${CLEAN_DIR}
for i in $(cat ./list); do
    fastp -i ${RAW_DIR}/${i}/${i}_1.fq.gz \
          -I ${RAW_DIR}/${i}/${i}_2.fq.gz \
          -w 16 \
          -o ${CLEAN_DIR}/${i}_R1.clean.fq.gz \
          -O ${CLEAN_DIR}/${i}_R2.clean.fq.gz \
          -h ${CLEAN_DIR}/${i}.html
done

# Post-trimming QC (commented out)
# mkdir clean-qc
# fastqc -o ./clean-qc ${CLEAN_DIR}/*

###================== 2. Alignment ===================
# Build index (commented out)
# bowtie2-build ${GENOME_REF} Anas_platyrhynchos.ASM874695v1.dna.toplevel.fa

# Alignment with bowtie2
mkdir -p samfile
for i in $(cat ./list); do
    bowtie2 -p 32 \
            --very-sensitive \
            -X 2000 \
            -x ${GENOME_REF} \
            -1 ${CLEAN_DIR}/${i}_R1.clean.fq.gz \
            -2 ${CLEAN_DIR}/${i}_R2.clean.fq.gz \
            -S ./samfile/${i}.sam
done

# SAM to BAM conversion and processing
mkdir -p bamfile sortbam
for i in $(cat ./list); do
    samtools view -@ 16 -q 30 -bF 4 -S ./samfile/${i}.sam -o ./bamfile/${i}.bam
    samtools sort -@ 16 -o ${BAM_DIR}/${i}.bam ./bamfile/${i}.bam
    samtools view -h ${BAM_DIR}/${i}.bam | grep -v chrM | samtools view -b -o ${BAM_DIR}/${i}_final.bam
done

###================== 3. Deduplication ===================
mkdir -p ${DEDUP_DIR}
for i in $(cat ./list); do
    picard MarkDuplicates \
        REMOVE_DUPLICATES=true \
        I=${BAM_DIR}/${i}_final.bam \
        O=${DEDUP_DIR}/${i}_deduplicate.bam \
        M=${DEDUP_DIR}/${i}.log
done

###================== 4. Peak Calling ===================
mkdir -p ${PEAK_DIR}
for i in $(cat ./list); do
    macs2 callpeak \
        -t ${DEDUP_DIR}/${i}_deduplicate.bam \
        -g 1211992756 \
        --nomodel \
        --shift -100 \
        --extsize 200 \
        -B \
        -n ${i} \
        --outdir ${PEAK_DIR}
done

###================== 5. Coverage Analysis ===================
# Index deduplicated BAMs
for i in $(cat ./list); do
    samtools index ${DEDUP_DIR}/${i}_deduplicate.bam
done

# Generate tracks
mkdir -p ${BW_DIR}
for i in $(cat ./list); do
    bamCoverage -bs 10 \
                --normalizeUsing RPGC \
                --effectiveGenomeSize 1211992756 \
                --extendReads \
                -b ${DEDUP_DIR}/${i}_deduplicate.bam \
                -o ${BW_DIR}/${i}.bw
done

###================== 6. TSS Analysis ===================
mkdir -p ${TSS_DIR}
for i in $(cat ./list); do
    computeMatrix reference-point \
        --referencePoint TSS \
        -p 15 \
        -b 3000 \
        -a 3000 \
        -R ${GENOME_BED} \
        -S ${BW_DIR}/${i}.bw \
        --skipZeros \
        -o ${TSS_DIR}/${i}_TSS.gz \
        --outFileSortedRegions ${TSS_DIR}/${i}.bed
done
