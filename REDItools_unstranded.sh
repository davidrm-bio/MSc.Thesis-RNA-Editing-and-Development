#!/bin/bash

# Parallel run # 
bin="/binf-isilon/rennie/gsn480/scratch/bin/reditools2.0-master"

source $bin/ENV/bin/activate


printf '\n Parallel Analysis on Samples'
date 


REF="/binf-isilon/rennie/gsn480/data/scratch/refGenome/MusMusculus/Mus_musculus.GRCm39.Pri.DNA.fa"
REF_SIZE="/binf-isilon/rennie/gsn480/data/scratch/refGenome/MusMusculus/Mus_musculus.GRCm39.Pri.DNA.fa.fai"
NUM_CORES=15
BAM_DIR="./"

for FILE in  "$BAM_DIR"/*.bam

do
	SAMPLE_NAME="$(basename "$FILE" .bam)"
	COV_FILE="${SAMPLE_NAME}.cov"
	COV_DIR="${SAMPLE_NAME}_coverage/"
	OUT_FILE="${SAMPLE_NAME}.out.gz"

	$bin/extract_coverage_dynamic.sh "$FILE" "$COV_DIR" "$REF_SIZE"
	mpirun -np "$NUM_CORES" $bin/src/cineca/parallel_reditools.py -f "$FILE" -o "$OUT_FILE" -S -s 0 -r "$REF" -bq 20 -q 20 -G "$COV_DIR$COV_FILE" -D "$COV_DIR" -t "${SAMPLE_NAME}_tmp/" -Z $REF_SIZE
	$bin/merge.sh "${SAMPLE_NAME}_tmp/" "$OUT_FILE" "$NUM_CORES"

done

echo '### Parallel Analysis finished ### '
deactivate
