#!/bin/bash

path="/storage-04/hic_GC/follicle_epi/dedup_bam/"

## binarized model
for i in `cat ./list`;
do
java -mx4000M -jar /storage-01/poultrylab1/yujz/software/chromHMM/ChromHMM/ChromHMM.jar BinarizeBam  -paired  ./cau_ASM874695v1.sizes ${path} ./${i}_marker_chromhmm.txt ${i}_chromhmm_binary
## learn model state 15
java -mx4000M -jar /storage-01/poultrylab1/yujz/software/chromHMM/ChromHMM/ChromHMM.jar LearnModel -holdcolumnorder -holdroworder ${i}_chromhmm_binary  ${i}_marker_chromhmm  15 mallard
done

