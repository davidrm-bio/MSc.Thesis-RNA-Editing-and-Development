#!/bin/bash

BAM_DIR="./"
THREADS=15
REGIONS="/binf-isilon/rennie/gsn480/data/scratch/refGenome/MusMusculus/mouse.bed"


for FILE in "$BAM_DIR"/*bam

do   
	SAMPLE_NAME="$(basename "$FILE" Aligned.sortedByCoord.out.bam)"
	FILE_TMP="${SAMPLE_NAME}tmp.bam"  
	FILE_DUPL="${SAMPLE_NAME}pre.rmdupl.bam"
	FINAL_FILE="${SAMPLE_NAME}clean.bam"
	
	# Remove unmapped reads; 
	# Keep only properly paired reads; 
	# Minimum mapping quality of 20;
	# Select only main chr (1-19 + X, Y)
	samtools view -@ "$THREADS" -F 0x04 -f 0x2 -q 20 -L $REGIONS -bh -o "$FILE_TMP" "$FILE"
	
	# Sort by query name ; Mark Duplicates
	samtools sort -@ "$THREADS" -n "$FILE_TMP" -T "$SAMPLE_NAME.tmp1" -o "$FILE_TMP"
	samtools fixmate -@ "$THREADS" -rm "$FILE_TMP" "$FILE_DUPL"

	# Sort by coordinate ; Remove Duplicates
	samtools sort -@ "$THREADS" -o "$FILE_DUPL" -T "$SAMPLE_NAME.tmp2" "$FILE_DUPL" 
	samtools markdup -@ "$THREADS" -r "$FILE_DUPL" "$FINAL_FILE"
	
	# Index files
	samtools index -@ "$THREADS" "$FINAL_FILE"   
    
	# Remove Intermediate
	rm "$FILE_TMP"
	rm "$FILE_DUPL"

done