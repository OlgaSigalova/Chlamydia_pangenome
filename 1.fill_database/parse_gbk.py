from Bio import SeqIO
import sys
import os
import pandas as pd
import numpy as np
from collections import Counter 
import re
import itertools
import fnmatch

################ Input: modify as needed ##########################################

gbk_dir = "/home/bsgroup/Chlamydia/genomes_apr_2019"
output_dir = "/home/bsgroup/Chlamydia/gbk_parsed"

# match GBK files in the input directory
gbk_files = fnmatch.filter(os.listdir(gbk_dir), '*.gbk')

###################################################################################


################ Output  ###########################################################

files=('Chlamydia_proteins.fasta','Chlamydia_proteins.txt', 'Chlamydia_frameshifts.txt', 'Chlamydia_frameshifted_proteins.txt')

for out_file in files:
    try:
        os.remove(output_dir + out_file)
    except:
        pass
        
output_handle1 = open(output_dir+'Chlamydia_proteins.fasta', "a")

output_handle2 = open(output_dir+'Chlamydia_proteins.txt', "a")
output_handle2.write("OrthMCL_prot_code\tOrthMCL_genome_code\tproduct\tstart\tend\tstrand\tlength\tframeshift\tnucl_seq\ttranslation\tGO_list\tEC_list\n")

output_handle3 = open(output_dir+'Chlamydia_frameshifts.txt', "a")
output_handle3.write("frameshift_id\tOrthMCL_prot_code\tfs_number\tto_end\tgap_size\n")

output_handle4 = open(output_dir+'Chlamydia_frameshifted_proteins.txt', "a")
output_handle4.write("OrthMCL_prot_code\tnum_frameshifts\tto_end_arr\tgap_size_arr\tnucl_seq_uc\n")

###################################################################################


r=0 # counter for genomes

for f in gbk_files:
    
    genome_id = f[:-4]
    input_handle  = open(gbk_dir + "/" + f, "r")
    
    gb_records = list(SeqIO.parse(input_handle, "genbank")) # relevant for gbk files with multiple contigs
    print str(r), genome_id
    
    r+=1
    k=0 # counter for CDSs per genome 
    n=0 # counter for frameshifts per genome
    
    for gb_record in gb_records:
        
        for seq_feature in gb_record.features :
            
            if seq_feature.type=="CDS" :
                
                k+=1
                
                # generate protein ID
                OrthMCL_prot_code = genome_id + "_" + str(k) 
                
                # check that translation is given on one line
                assert len(seq_feature.qualifiers['translation']) == 1 
                
                # Basic info
                
                start = seq_feature.location.start.position
                end = seq_feature.location.end.position
                strand = seq_feature.location.strand
                product = seq_feature.qualifiers['product'][0] # annotation
                nucl_seq = seq_feature.extract(gb_record.seq) # sequences after correction for frameshifts
                aa_seq=seq_feature.qualifiers['translation'][0] 
                lenght = len(aa_seq)

                # Frameshifts
                
                frameshift = "no"
                
                if seq_feature.location_operator == 'join':
                    
                    frameshift = "yes"
                    num_frameshifts = len(seq_feature.sub_features) - 1
                    gap_size_arr = [] # distance between joined fragments. if negative, the fragment overlap
                    to_end_arr = [] # number of nucleotides from frameshift to gene end
                    
                    if strand == -1:
                        
                        nucl_seq_uc = gb_record.seq[start:end].reverse_complement() # uncorrected nc sequences
                        
                        for i in range(num_frameshifts):
                            
                            n+=1
                            frameshift_id = genome_id + "_fs"+str(n)
                            fs_number = num_frameshifts-i
                            gap_size = str(seq_feature.sub_features[i+1].location.start-seq_feature.sub_features[i].location.end)
                            gap_size_arr.append(gap_size)
                            to_end=str(seq_feature.sub_features[i].location.end-seq_feature.sub_features[0].location.start)
                            to_end_arr.append(to_end) 
                            
                            # write summary per frameshift
                            line = "\t".join([frameshift_id, OrthMCL_prot_code, str(fs_number), str(to_end) , str(gap_size)])
                            output_handle3.write(line + "\n")
                            
                    else:
                        
                        nucl_seq_uc=gb_record.seq[start:end] # uncorrected nc sequences
                        
                        for i in range(num_frameshifts):
                            n+=1
                            frameshift_id = genome_id+"_fs"+str(n)
                            fs_number = i+1
                            gap_size = str(seq_feature.sub_features[i+1].location.start-seq_feature.sub_features[i].location.end)
                            gap_size_arr.append(gap_size)
                            to_end = str(seq_feature.sub_features[num_frameshifts].location.end-seq_feature.sub_features[i+1].location.start)
                            to_end_arr.append(to_end)
                            
                            # write summary per frameshift
                            line = "\t".join([frameshift_id, OrthMCL_prot_code, str(fs_number), str(to_end) , str(gap_size)])
                            output_handle3.write(line + "\n")
                            
                    gap_size_arr=";".join(gap_size_arr)
                    to_end_arr=";".join(to_end_arr)
                    
                    # Write summary of frameshifts per protein
                    line = "\t".join([OrthMCL_prot_code, str(num_frameshifts) , str(to_end_arr), str(gap_size_arr), str(nucl_seq_uc)])
                    output_handle4.write(line + "\n")

                    
                #get EC numbers 
                
                EC=""
                
                try:
                    EC = seq_feature.qualifiers["EC_number"]
                    EC = ";".join(EC)
                except:
                    pass
                #get GOs
                GOs = ""
                if len(seq_feature.qualifiers['db_xref']) > 1:
                    GOs = seq_feature.qualifiers['db_xref'][1:]
                    GOs = ";".join(GOs)
                    

                # write protein fasta
                output_handle1.write(">%s\n%s\n" % (OrthMCL_prot_code, str(aa_seq)))
                
                # write protein info
                line = "\t".join([OrthMCL_prot_code, genome_id, product, str(start), str(end), str(strand), str(lenght), frameshift, str(nucl_seq), str(aa_seq),  GOs, EC])
                output_handle2.write(line + "\n")
                
                


output_handle1.close()
output_handle2.close()
output_handle3.close()
output_handle4.close()