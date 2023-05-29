# MSc.Thesis-RNA-Editing-and-Development

## Part 1 - Benchmarking of Bioinformatic tools
Several RNA editing detection tools were used to analysed RNA-seq data to assess their performance. The bioinformatic tools included are: 1) RED-ML, 2) SPRINT, 3) REDItools2, 4) JACUSA2 and 5) BCFtools. Additionally, three different aligners are used (BWA, HISAT2 and STAR). 

Scripts can be found in folder "Benchmarking". Available scripts:

-  Pre-processing: Bash scripts used to pre-processed FASTQ files, map to the reference genome and do quality control on BAM files.
- Tools: Bash scripts used to run the analysis with the bioinformatic tools
-  Downstream: Python and R scripts used in the downstream processing of the results
- Databases: Python scripts used to convert to JSON format files from dbSNP (VCF file) and REDIportal database.


## RNA editing in the development of cardiomyocytes
REDItools2 was used to analysed RNA-seq data to study the involvement of RNA editing in the development of cardiomyocytes

- Mapping: Bash script used to pre-processed FASTQ files and map to the reference Genome using STAR.
- PostMapping: Bash scripts used to process BAM files.
- REDItools2: Bash script used to run REDItools2.
- Downstream: collection of Python and R scripts used in the downstream analysis


