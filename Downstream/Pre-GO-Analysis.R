# Libraries
library (GenomicRanges)
library (tidyverse)
library(rtracklayer)

# Part 1 - Loading data
gtf <- import('/Users/david/Downloads/DATA/MusMusculus.GRCm39.Pri.OnlyGenesAnnotation.gtf') 
bed <- import ('MiceEP-WT>KO.bed') ## Change accordingly
bed2 <- import ('MiceEP-WT<KO.bed') ## Change accordingly

# Part 2 - Analysis 
results <- as.data.frame (findOverlaps (bed,gtf))
results2 <- as.data.frame (findOverlaps (bed2,gtf))


# Part 3 - Processing
clean_results <- function (df, result){
  df %>% 
    as.data.frame() %>% 
    rownames_to_column('Count') %>% 
    filter(Count %in% result$subjectHits) %>% 
    filter (gene_type == 'protein_coding') %>% 
    select(seqnames, strand, gene_name) %>% 
    select(gene_name) -> df.clean
  return(df.clean)
}

results.clean <- clean_results(gtf, results)
results2.clean <- clean_results(gtf, results2)


# Part 4 - Saving
write.table (results.clean, 'MiceEP-WT>KO.txt',  quote = FALSE, row.names = F, col.names = F)
write.table (results2.clean, 'MiceEP-WT<KO.txt',  quote = FALSE, row.names = F, col.names = F)
