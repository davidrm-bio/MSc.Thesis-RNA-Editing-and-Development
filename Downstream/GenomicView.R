# Library
library (GenomicRanges)
library (tidyverse)
library(rtracklayer)


# Part 1 - Loading data
gtf <- import('/Users/david/Downloads/DATA/MusMusculus.GRCm39.Pri.OnlyGenesAnnotation.gtf') 

srr11 <- import ('Mice_KO/SRR10581811_clean.bed')
srr12 <- import ('Mice_KO/SRR10581812_clean.bed')
srr13 <- import ('Mice_KO/SRR10581813_clean.bed')
srr14 <- import ('Mice_KO/SRR10581814_clean.bed')


# Part 2 - Select a gene
gtf %>% 
  as.data.frame() %>% 
  filter(gene_name == 'Flna') %>%  # Change accordingly
  makeGRangesFromDataFrame ()  -> gene


# Part 3 - Analysis
out11 <- as.data.frame (findOverlaps (srr11,gene))
out12<- as.data.frame (findOverlaps (srr12,gene))
out13<- as.data.frame (findOverlaps (srr13,gene))
out14<- as.data.frame (findOverlaps (srr14,gene))

clean_data <- function(df) {
    df %>% 
  as.data.frame() %>% 
  rownames_to_column('Count') %>% 
  filter(Count %in% out11$queryHits) -> df.clean

  return (df.clean)
}

save11 <- clean_data(srr11)
save12 <- clean_data(srr12)
save13 <- clean_data(srr13)
save14 <- clean_data(srr14)

# Part 4 - Save to file
write2Bed <- function (df){
  df %>% 
    mutate (start = end-1) %>% 
    mutate (end = end) %>% 
    select(seqnames, start, end, name, score, strand) -> df.clean
  return (df.clean)
}

write.table(write2Bed(save11), 'USCS-SRR11-genes.bed', quote = F, col.names = F, row.names = F, append=T)
write.table(write2Bed(save12), 'USCS-SRR12-gnes.bed', quote = F, col.names = F, row.names = F, append=T)

write.table(write2Bed(save13), 'USCS-SRR13-genes.bed', quote = F, col.names = F, row.names = F, append=T)
write.table(write2Bed(save14), 'USCS-SRR14-genes.bed', quote = F, col.names = F, row.names = F, append=T)
