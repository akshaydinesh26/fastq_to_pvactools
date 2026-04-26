setwd("D:/Users/Akshay.Dinesh/Desktop/nextflow_workflow/viral_workflow/lut/")

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

data <- fread("mart_export_human.txt")
head(data)
names(data) <- c("gene_id_strip","gene_id","transcript_id_strip",
                 "transcript_id","protein_id","chr","start","end","strand",
                 "transcript_length","gene_name","gene_synonym",
                 "gene_type")

final <- data %>% 
  mutate(transcript_length=transcript_length/1000) %>%
  dplyr::select(gene_id,gene_name,gene_type,
         strand,chr,start,end,transcript_id,transcript_length,gene_id)

final <- final %>%
  mutate(original_gene_id=gene_id) %>%
  group_by(gene_id,gene_name,gene_type,strand,original_gene_id,chr,start,
           end) %>%
  summarise(median_transcript_length_kb=median(transcript_length),
            transcript_id=paste0(transcript_id,collapse = ","),
            transcript_length=paste0(transcript_length,collapse = ",")) %>%
  ungroup()

final <- final %>% dplyr::select(gene_id,gene_name,gene_type,median_transcript_length_kb,strand,chr,start,end,
                          transcript_id,transcript_length,original_gene_id)

final <- final %>% mutate(strand=ifelse(strand==1,"+",ifelse(strand==-1,"-",strand)))

head(final)


chr_valid <- c(1:22,"X","Y","MT")
final_chr <- final %>% filter(chr %in% chr_valid)
unique(final_chr$chr)
ids_in_ref <- fread("SRR28960393_ReadsPerGene.out.tab",skip = 4) 
names(ids_in_ref) <- c("gene_id_old","unstraanded","stranded","reverse")


ids_in_ref <- ids_in_ref %>% mutate(index=str_remove(gene_id_old,"\\..*"))
final_chr <- final_chr %>% mutate(index=str_remove(gene_id,"\\..*"))
d <- final_chr %>% left_join(ids_in_ref,by="index") %>% select(index,gene_id_old,gene_id) %>%
View(final_chr)
write.table(final_chr,"reference.LUT.human.gtf",sep="\t",quote=T,row.names = F)
