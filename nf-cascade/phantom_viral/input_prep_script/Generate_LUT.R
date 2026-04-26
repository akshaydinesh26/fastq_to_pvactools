setwd("D:/Users/Akshay.Dinesh/Desktop/nextflow_workflow/viral_workflow/lut/")


library(data.table)
library(dplyr)
library(tidyr)
library(stringr)


data <- fread("data_table.txt")
head(data)
names(data) <- c("chr","start","end","strand","gene_id","gene_name","gene_type",
                 "transcript_id","protein_id")
data <- data %>% mutate(gene_id=str_split_i(gene_id,'"',2),
                        gene_name=str_split_i(gene_name,'"',2),
                        gene_type=str_split_i(gene_type,'"',2),
                        transcript_id=str_split_i(transcript_id,'"',2),
                        protein_id=str_split_i(protein_id,'"',2))


final <- data %>% filter(gene_type=="protein_coding") %>% 
  mutate(transcript_length=(end-start)/1000) %>%
  select(gene_id,gene_name,gene_type,
         strand,chr,start,end,transcript_id,transcript_length,protein_id,gene_id)

final <- final %>%
  mutate(virus=chr) %>%
  group_by(gene_id,gene_name,gene_type,strand,virus,chr,start,
           end) %>%
  summarise(median_transcript_length_kb=median(transcript_length),
            transcript_id=paste0(transcript_id,collapse = ","),
            transcript_length=paste0(transcript_length,collapse = ","),
            protein_id=paste0(protein_id,collapse = ",")) %>%
  ungroup()


final <- final %>% select(gene_id,gene_name,gene_type,median_transcript_length_kb,strand,chr,start,end,
                          transcript_id,transcript_length,protein_id,virus)

final <- final %>% distinct()
final <- final %>% filter(protein_id != "")
final %>% group_by(gene_id) %>% filter(n()>1) %>% ungroup() %>% View()
#length(unique(final$gene_id))


write.table(final,"reference.LUT.virus.gtf",sep="\t",quote=F,row.names = F)
