#!/bin/bash

# Activate Conda environment (commented out)
conda activate hicpro3.1

###================== 1. HiC-Pro Analysis ===================
# Note: HiC-Pro pipeline command (commented out for reference)
 /storage-01/poultrylab1/lizhen/pk_hic_shice/hicpro/HiC-Pro/bin/HiC-Pro \
   -i ./raw_data/ \
   -o result \
   -c config-hicpro.txt

###================== 2. Compartment Identification ===================
# Define reference genome path
DUCK_GENOME=/storage-02/lizhen/genome/Anas_platyrhynchos.ASM874695v1.dna.toplevel.fa
HOMER_DIR=./homer

## 2.1 Create Tag Directory
# Note: Command to create HOMER tag directory (commented out)
 makeTagDirectory ${HOMER_DIR}/SYF_homer \
   -format HiCsummary \
   SYF.allValidPairs.homer

## 2.2 Calculate Principal Components
# Note: Commands to run Hi-C PCA (commented out)
 runHiCpca.pl 40k ${HOMER_DIR}/SYF_homer \
   -cpu 16 \
   -res 40000 \
   -genome ${DUCK_GENOME} \
   -pc 1

# Output: Generates bedGraph and txt files for IGV visualization

###================== 3. TAD Identification ===================
## 3.1 Convert HiC-Pro matrix to h5 format
# Note: Conversion command (commented out)
hicConvertFormat \
   -m /storage-01/poultrylab1/lizhen/pk_hic_shice/f5_result/hic_results/matrix/F5/iced/20000/F5_20000_iced.matrix \
   --inputFormat hicpro \
   -bf /storage-01/poultrylab1/lizhen/pk_hic_shice/f5_result/hic_results/matrix/F5/raw/20000/F5_20000_abs.bed \
   --outputFormat h5 \
   --outFileName F5_20K.h5

## 3.2 Identify TADs
# Note: TAD calling command (commented out)
 hicFindTADs \
   -m F5_20K.h5 \
   --outPrefix F5_20K_thres0.05_delta0.01_fdr \
   --thresholdComparisons 0.05 \
   --delta 0.01 \
   --correctForMultipleTesting fdr \
   -p 30
