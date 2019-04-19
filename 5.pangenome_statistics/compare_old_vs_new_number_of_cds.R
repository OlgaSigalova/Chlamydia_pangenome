
# old CDS stat
path_cds = "/Users/sigalova/Desktop/Chlamydia_pangenome/info/genome_n_cds_by_contig_raw_apr2019.csv"
df_cds = read.csv(path_cds)

# contigs info
path_contigs = "/Users/sigalova/Desktop/Chlamydia_pangenome/info/genome_info_by_contig_final_apr2019.csv"
df_contigs = read.csv(path_contigs)

df = merge(df_cds, df_contigs, by = "accession")

# summarize by gene
df %<>% group_by(genome_id) %>% summarize(n_cds_old = sum(n_cds))

# new CDS stat
path = "/Users/sigalova/Desktop/Chlamydia_pangenome/data/gbk_parsed/filt/Chlamydia_proteins.txt"
proteins = fread(path, sep = "\t")

proteins_sum = proteins %>% 
  group_by(OrthMCL_genome_code) %>% 
  summarize(n_cds_new = n()) %>% 
  rename(genome_id = OrthMCL_genome_code)

df = merge(df, proteins_sum, by = "genome_id") 
df %<>% arrange(-n_cds_old)
 
# on average, 40 extra CDS in new annotation
df %>% filter(n_cds_old > 0) %>% summarize(mean(n_cds_new - n_cds_old))

outf = "/Users/sigalova/Desktop/Chlamydia_pangenome/info/old_vs_new_number_of_cds.csv"
write.csv(df, outf, row.names = FALSE)
