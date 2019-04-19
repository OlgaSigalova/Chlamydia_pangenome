library(tidyverse)
library(magrittr)
options(stringsAsFactors = FALSE)

###### finalize and merge with old table - after manual correction of species/strains ######

f_new = "/Users/sigalova/Desktop/Chlamydia_pangenome/info/genome_info_final_sp_corr_apr2019.csv"
df_new = read.csv(f_new) %>% 
  select(genome_id, accession, full_name, species, strain, genome_length, n_contigs)

f_old = "/Users/sigalova/Desktop/Chlamydia_pangenome/info/genomes_info_jun_2016.csv"
used_genomes = "/Users/sigalova/Desktop/Chlamydia_pangenome/info/final_list_of_genomes_jun2016.csv"
df_old = read.csv(f_old)
gl = read.csv(used_genomes)$genome_id
df_old %<>% filter(genome_id %in% gl) %>%
  rename(accession = ref, n_contigs = contigs, genome_length = DNAlength) %>%
  select(genome_id, accession, full_name, species, strain, genome_length, n_contigs)


df_final = rbind.data.frame(df_old, df_new) %>% arrange(genome_id)

outf = "/Users/sigalova/Desktop/Chlamydia_pangenome/info/genome_info_combined.csv"
write.csv(df_final, outf, row.names = FALSE)
