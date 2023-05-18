#!/bin/bash

NCUT=11
BASE='scratch' 
TMP='tmp.dir'
PROCESS='fastq_process'
THR=10 ## threads


date
printf '\nSTAR RNA-seq mapping########## \n'
STAR --version


### STAR Genome Index

GENDIR="/binf-isilon/rennie/gsn480/data/scratch/refGenome/MusMusculus/STAR_135"

FASTQS=$( find ./scratch/fastq/ -type f -name "*fastq.gz" -exec ls {} + )
SAMPVEC=$( ls $FASTQS | sed 's!.*/!!' | rev | cut -c ${NCUT}- | rev  | sort | uniq )
SAMPNUM=$( ls ${FASTQS} | sed 's!.*/!!' | rev | cut -c ${NCUT}- | rev  | sort | uniq | wc -l)

printf "SAMPLES DETECTED = ${SAMPNUM} \n"
printf "FASTQ'S DETECTED = $( ls $FASTQS | wc -l ) \n"


for SAMP in $( echo $SAMPVEC );

do  
    printf "\n"
    echo "SAMPLE##################"
    
    echo $SAMP
    R1=$( ls $FASTQS| grep ${SAMP} | grep "_1" )
    R2=$( ls $FASTQS| grep ${SAMP} | grep "_2" )
    echo "READ1####"
    ls -lh $R1
    echo "READ2####"
    ls -lh $R2

    MR1="${BASE}/${TMP}/Read1_${SAMP}.fastq"
    MR2="${BASE}/${TMP}/Read2_${SAMP}.fastq"

    cat $R1 | gunzip > $MR1
    cat $R2 | gunzip > $MR2

    echo "MERGED ####"
    ls -lh $MR1
    ls -lh $MR2

    TR1="${BASE}/${PROCESS}/Read1_Processed_${SAMP}.fastq"
    TR2="${BASE}/${PROCESS}/Read2_Processed_${SAMP}.fastq"

    fastp -i $MR1 -I $MR2 -o $TR1 -O $TR2 --cut_front --cut_tail --detect_adapter_for_pe  --overrepresentation_analysis --correction --trim_front1 7 --trim_front2 7 -h ${SAMP}report.html -w ${THR}

    echo "PROCCESSED ##########"
    ls -lh $TR1
    ls -lh $TR2

    echo 'STAR'
    STAR --runMode alignReads --genomeDir $GENDIR --readFilesIn $TR1 $TR2 --runThreadN $THR --outSAMstrandField intronMotif --outSAMmode Full --outSAMattributes All --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix ./${BASE}/mapping/${SAMP} --outTmpDir ./${BASE}/mapping/tmp${SAMP}
    
done

echo 'Job Finished'