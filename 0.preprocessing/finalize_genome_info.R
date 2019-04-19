library(tidyverse)
library(magrittr)
options(stringsAsFactors = FALSE)

# genomes before RAST - by contig
f_raw = "/Users/sigalova/Desktop/Chlamydia_pangenome/info/genome_info_by_contig_raw_apr2019.csv"
# genomes after RAST and filtering - by contig
f_annot = "/Users/sigalova/Desktop/Chlamydia_pangenome/info/genome_info_by_contig_annot_filt_apr2019.csv"

df_raw = read.csv(f_raw)
df_annot = read.csv(f_annot)

df = merge(df_raw, df_annot, by = "accession")

# genomes final - by contig
outf_by_contif = "/Users/sigalova/Desktop/Chlamydia_pangenome/info/genome_info_by_contig_final_apr2019.csv"
write.csv(df, outf_by_contif, row.names = FALSE)

# genomes final - by genome
outf_by_genome= "/Users/sigalova/Desktop/Chlamydia_pangenome/info/genome_info_final_apr2019.csv"

df_sum = df %>% 
  group_by(genome_id) %>% 
  mutate(max_length = max(DNA_length), 
         genome_length = sum(DNA_length),
         n_contigs = n()) %>%
  filter(DNA_length == max_length) %>%
  ungroup() %>%
  select(-DNA_length, -max_length) 

# draft species and strain to be corrected
df_sum$species = unlist(lapply(df_sum$full_name, function(x) strsplit(x, " ")[[1]][2]))
df_sum$strain = unlist(lapply(df_sum$full_name, function(x) paste(strsplit(x, " ")[[1]][3:5], collapse = " ")))

df_sum$strain = gsub("strain ", "", df_sum$strain)
df_sum$strain = gsub("isolate ", "", df_sum$strain)
df_sum$strain = gsub(" genome", "", df_sum$strain)
df_sum$strain = gsub(" isolate", "", df_sum$strain)
df_sum$strain = gsub(" chromosome.*", "", df_sum$strain)
df_sum$strain = gsub(".contig.*", "", df_sum$strain)
df_sum$strain = gsub(", complete", "", df_sum$strain)
df_sum$strain = gsub(" sequence.", "", df_sum$strain)
df_sum$strain = gsub(". NA", "", df_sum$strain)
df_sum$strain = gsub(" assembly.*", "", df_sum$strain)
df_sum$strain = gsub("corallus ", "", df_sum$strain)
df_sum$strain = gsub(" ", "_", df_sum$strain)


write.csv(df_sum, outf_by_genome, row.names = FALSE)
