#!/bin/bash

# Define directory paths
RAW_DIR=/storage-02/duck_library/follicle/cut_tag
CLEAN_DIR=./cleandata
SAM_DIR=./samfile
BAM_DIR=./bamfile
SORT_BAM_DIR=./sortbam
DEDUP_DIR=./dedup_bam
BW_DIR=./bw
TSS_DIR=./tss
PEAK_DIR=./callpeak
BOWTIE_SUMMARY=./bowtie2_summary
HTML_DIR=./html
GENOME_REF=./Anas_platyrhynchos.ASM874695v1.dna.toplevel.fa
GENOME_BED=./Anas_platyrhynchos.ASM874695v1.110.bed
EFFECTIVE_GENOME_SIZE=1211992756

###================== 1. Quality Control ===================
mkdir -p ${CLEAN_DIR} ${HTML_DIR}
for i in $(cat ./list); do
    fastp -i ${RAW_DIR}/${i}/${i}_1.fq.gz \
          -I ${RAW_DIR}/${i}/${i}_2.fq.gz \
          -w 32 \
          -o ${CLEAN_DIR}/${i}_R1.clean.fq.gz \
          -O ${CLEAN_DIR}/${i}_R2.clean.fq.gz \
          -h ${HTML_DIR}/${i}.html
done

###================== 2. Alignment ===================
mkdir -p ${SAM_DIR} ${BOWTIE_SUMMARY} ${BAM_DIR} ${SORT_BAM_DIR}
for i in $(cat ./list); do
    bowtie2 --end-to-end \
            --very-sensitive \
            --no-mixed \
            --no-discordant \
            -p 32 \
            -I 10 \
            -X 700 \
            -x ${GENOME_REF} \
            -1 ${CLEAN_DIR}/${i}_R1.clean.fq.gz \
            -2 ${CLEAN_DIR}/${i}_R2.clean.fq.gz \
            -S ${SAM_DIR}/${i}.sam \
            > ${BOWTIE_SUMMARY}/${i}_bowtie2.txt
done

# Convert SAM to BAM and sort
for i in $(cat ./list); do
    samtools view -@ 32 -q 25 -bF 4 -S ${SAM_DIR}/${i}.sam -o ${BAM_DIR}/${i}.bam
    samtools sort -@ 32 -o ${SORT_BAM_DIR}/${i}.bam ${BAM_DIR}/${i}.bam
done

###================== 3. Duplicate Processing ===================
mkdir -p ${DEDUP_DIR}
for i in $(cat ./list); do
    # Mark duplicates without removal
    picard MarkDuplicates \
        I=${SORT_BAM_DIR}/${i}.bam \
        O=${DEDUP_DIR}/${i}_dupMarked.bam \
        M=${DEDUP_DIR}/${i}_dupMarked.txt
    
    # Remove duplicates
    picard MarkDuplicates \
        REMOVE_DUPLICATES=true \
        I=${SORT_BAM_DIR}/${i}.bam \
        O=${DEDUP_DIR}/${i}_deduplicate.bam \
        M=${DEDUP_DIR}/${i}_rmdep.txt
done

###================== 4. Coverage Analysis ===================
# Index deduplicated BAMs
for i in $(cat ./list); do
    samtools index ${DEDUP_DIR}/${i}_deduplicate.bam
done

# Generate coverage tracks
mkdir -p ${BW_DIR}
for i in $(cat ./list); do
    bamCoverage -bs 10 \
                --normalizeUsing RPGC \
                --effectiveGenomeSize ${EFFECTIVE_GENOME_SIZE} \
                --extendReads \
                -b ${DEDUP_DIR}/${i}_deduplicate.bam \
                -o ${BW_DIR}/${i}.bw
done

###================== 5. TSS Analysis ===================
mkdir -p ${TSS_DIR}
for i in $(cat ./list); do
    computeMatrix reference-point \
        --referencePoint TSS \
        -p 32 \
        -b 3000 \
        -a 3000 \
        -R ${GENOME_BED} \
        -S ${BW_DIR}/${i}.bw \
        --skipZeros \
        -o ${TSS_DIR}/${i}_TSS-3k.gz \
        --outFileSortedRegions ${TSS_DIR}/${i}.bed
done

# Generate visualization plots
for i in $(cat ./list); do
    plotHeatmap -m ${TSS_DIR}/${i}_TSS-3k.gz \
                -out ${TSS_DIR}/${i}_Heatmap-3k.pdf \
                --plotFileFormat pdf \
                --dpi 720
    
    plotProfile -m ${TSS_DIR}/${i}_TSS-3k.gz \
                -out ${TSS_DIR}/${i}_Profile-3k.png \
                --dpi 720
done

###================== 6. Peak Calling ===================
mkdir -p ${PEAK_DIR}

# Broad peaks (H3K4me1, H3K27me3)
for i in $(cat ./kuan.txt); do
    macs2 callpeak \
        -t ${DEDUP_DIR}/${i}_deduplicate.bam \
        --broad \
        --broad-cutoff 0.1 \
        -g ${EFFECTIVE_GENOME_SIZE} \
        -f BAMPE \
        -q 0.1 \
        -n ${i} \
        --outdir ${PEAK_DIR}
done

# Narrow peaks (H3K4me3, H3K27ac)
for i in $(cat ./zhai.txt); do
    macs2 callpeak \
        -t ${DEDUP_DIR}/${i}_deduplicate.bam \
        -g ${EFFECTIVE_GENOME_SIZE} \
        -q 0.1 \
        -B \
        -n ${i} \
        --outdir ${PEAK_DIR}
done
