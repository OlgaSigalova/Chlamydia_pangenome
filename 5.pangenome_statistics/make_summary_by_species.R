library(tidyverse)
library(magrittr)
library(data.table)
options(stringsAsFactors = FALSE)


path = "/Users/sigalova/Desktop/Chlamydia_pangenome/info/genome_info_combined.csv"
df = read.csv(path)

path = "/Users/sigalova/Desktop/Chlamydia_pangenome/data/gbk_parsed/filt/Chlamydia_proteins.txt"
proteins = fread(path, sep = "\t")

proteins_sum = proteins %>% 
  group_by(OrthMCL_genome_code) %>% summarize(n_cds = n(), n_frameshifts = sum(frameshift == "yes")) %>% 
  rename(genome_id = OrthMCL_genome_code)

df = merge(df, proteins_sum, by = "genome_id")


df_sum = df %>% group_by(species) %>% 
  summarize(n_genomes = n(), 
            med_size = round(median(genome_length)), 
            med_n_cds = round(median(n_cds))) %>%
  arrange(-n_genomes)


outf = "/Users/sigalova/Desktop/Chlamydia_pangenome/info/dataset_summary.csv"
write.csv(df_sum, outf, row.names = FALSE)
