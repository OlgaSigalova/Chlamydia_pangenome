# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 11:59:21 2015

@author: olga
"""

from Bio import SeqIO
import sys
import os
#os.chdir("/home/olga/Public/Pangenome project/working scripts/")



# k - counters for proteins within genome

input_dir = sys.argv[1]
output_dir=sys.argv[2]
info_dir=sys.argv[3]

#info_dir="/home/olga/Public/Pangenome_project/"
# read taxonomy data
spp=dict()
stains=dict()
output_handle0 = open(info_dir+'genomes_info.csv', "r")
for line in output_handle0:
    line=line.rstrip().split(",")
    gi=line[1]
    stain=line[2]
    species=" ".join(stain.split(" ")[0:2])
    spp[gi]=species
    stains[gi]=stain
output_handle0.close() 




gbk_filenames=os.listdir(input_dir)

files=('Chlamydia_genomes.txt','Chlamydia_proteins.fasta','Chlamydia_proteins.txt','Chlamydia_protein_GOs.txt')

for out_file in files:
    try:
        os.remove(output_dir+out_file)
    except:
        pass
    
output_handle1 = open(output_dir+'Chlamydia_genomes.txt', "a")
output_handle1.write("OrthMCL_genome_code\tgi\tref\tspecies\tstain\tDNAlength\tDNA\n")
output_handle2=open(output_dir+'Chlamydia_proteins.fasta', "a")
output_handle3=open(output_dir+'Chlamydia_proteins.txt', "a")
output_handle3.write("OrthMCL_prot_code\tOrthMCL_genome_code\tproduct\tstart\tend\tstrand\tlength\tframeshift\tnucl_seq\tnucl_seq_corrected\ttranslation\n")
output_handle4=open(output_dir+'Chlamydia_protein_GOs.txt', "a")

#gbk_filename='/home/olga/Public/Pangenome_project/genomes_gbk/C113.gbk'

for gbk_filename in gbk_filenames:
    
    # read gbk file
    OrthMCL_genome_code=gbk_filename[0:4]
    input_handle  = open(input_dir+gbk_filename, "r") # read each gbk file 
    gb_record = SeqIO.read(input_handle, "genbank")
    print "Dealing with GenBank record %s" % gb_record.name
    
    #parse genome information
    info=gb_record.name.split("|")
    DNAlength=len(gb_record.seq) 
    assert info[0]=="gi" and info[2]=="ref" # check that gi and ref codes are provided
    gi=info[1]
    species=spp[gi]
    stain=stains[gi]
    ref=info[3]
    genome=gb_record.seq
    k=0
    line=str(OrthMCL_genome_code)+"\t"+str(gi)+"\t"+str(ref)+"\t"+species+"\t"+stain+"\t"+str(DNAlength)+"\t"+str(genome)
    output_handle1.write(line+"\n")
    
    # parse protein information
    
#    for gb_record in SeqIO.parse(input_handle, "genbank") :
       
    for seq_feature in gb_record.features :
        
        if seq_feature.type=="CDS" :
            k+=1
            OrthMCL_prot_code=OrthMCL_genome_code+"_"+str(k)
            OrthMCL_code=OrthMCL_genome_code+"|"+OrthMCL_prot_code
            assert len(seq_feature.qualifiers['translation'])==1 # check that translation is given on one line
            translation=seq_feature.qualifiers['translation'][0]            
            lenght=len(translation)
            # write protein sequences to file
            output_handle2.write(">%s\n%s\n" % (
                OrthMCL_code,
                translation))
                
            #check for joins
            if seq_feature.location_operator=='join':
                frameshift="yes"
            else:
                frameshift="no"
            
            #get EC numbers 
            EC="-"
            try:
#                seq_feature.qualifiers["EC_number"]
                EC=seq_feature.qualifiers["EC_number"]
            except:
                pass
            
            #get GOs

            if len(seq_feature.qualifiers['db_xref'])>1:
                GOs=seq_feature.qualifiers['db_xref'][1:]
                GOs=[x[3:] for x in GOs]
                GOs=";".join(GOs)
                output_handle4.write(" ".join([OrthMCL_code, GOs,"\n"]))  
            
            start=seq_feature.location.start.position
            end=seq_feature.location.end.position
            strand=seq_feature.location.strand
            product=seq_feature.qualifiers['product'][0]
            
            # the code below ignores joins
            if strand==-1:
                nucl_seq=gb_record.seq[start:end].reverse_complement()
            else:
                nucl_seq=gb_record.seq[start:end]
            
            # here the frameshifts are fixed
            nucl_seq_corrected=seq_feature.extract(gb_record.seq)
            
            output_handle3.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                OrthMCL_prot_code,
                OrthMCL_genome_code,
                product,
                start,
                end,
                strand,
                lenght,
                frameshift,
                nucl_seq,
                nucl_seq_corrected,
                translation))
                

output_handle1.close()
output_handle2.close()
output_handle3.close()
output_handle4.close()
print "Done"
