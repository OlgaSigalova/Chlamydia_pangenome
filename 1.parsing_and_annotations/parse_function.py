# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 18:40:35 2015

@author: olga
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 14:36:17 2015

@author: olga
"""
import os
from collections import Counter
import sys 
import re

in1 = sys.argv[1]
in2=sys.argv[2]
in3=sys.argv[3]
in4=sys.argv[4]
out1=sys.argv[5]

in_f1 = open(in1, "r")
in_f2 = open(in2, "r")
in_f3 = open(in3, "r")
in_f4 = open(in4, "r")
tmp=in_f1.readline()
tmp=in_f2.readline()
tmp=in_f3.readline()

#in1='/home/olga/Public/work/Pangenome_project/out/'
#in_f1 = open(in1+'gene_products.txt', "r")
#in2='/home/olga/Public/work/Pangenome_project/out/'
#in_f2 = open(in2+'OGgroups.txt', "r")
#out1='/home/olga/Public/work/Pangenome_project/out/'

try:
    os.remove(out1+'products_by_OGgroup.txt')
except:
    pass
out_f1 = open(out1+'products_by_OGgroup.txt', "a")

prod_dict=dict()
frame_dict=dict()
genome_dict=dict()

# proteins info
for line in in_f1:
    line=line.rstrip().split("\t")
    [a,b,c,d]=[line[0],line[2],line[7],line[1]]
    prod_dict[a]=b
    frame_dict[a]=c
    genome_dict[a]=d

species_dict=dict()

#genomes info    
for line in in_f2:
    line=line.rstrip().split("\t")
    [a,b]=[line[0],line[3].split(" ")[1]]
    species_dict[a]=b

species=set(species_dict.values())

genomes=species_dict.keys()
genomes.sort()

#GO info

GOdict=dict()

for line in in_f3:
    [a,b]=line.rstrip().split(" ")
    GOdict[a]=b.strip()

line="group_id\tgroup_size\tdistribution\tOGclass\tall_hypothetical\tproduct\twith_frameshift\tproduct_distr\ttrachomatis\tpsittaci\tpneumoniae\tpecorum\tabortus\tcaviae\tfelis\tmuridarum\t"+"\t".join(genomes)+"\tGO_distr\tGOset"+"\n"    
out_f1.write(line)

proteins_to_OG=dict()
    
#OG info
for line in in_f4:
    [group,proteins]=line.rstrip().split(":")

    proteins=proteins.strip().split(" ")
    OGsize=len(proteins)
    prod_list=[]
    frame_list=[]
    genome_list=[]
    GOlist=[]
    trachomatis, psittaci, pneumoniae, pecorum, abortus, caviae, felis, muridarum=[0,0,0,0,0,0,0,0]
    
    product=''
    distribution=0
    counts=[0]*len(genomes)
    OGclass='non_universal'
    for protein in proteins:
        
        try: 
            GOs=GOdict[protein].split(";")
            for GO in GOs:
                GO="GO:"+GO
                GOlist.append(GO)
        except:
            pass        
        
        protein=protein.split("|")[1]
        proteins_to_OG[protein]=group
        prod=prod_dict[protein]
        genome=genome_dict[protein]
        counts[int(genome[1:])-1]+=1
        species=species_dict[genome]
        exec(species + "+= 1")
        frame=frame_dict[protein]
        prod_list.append(prod)
        frame_list.append(frame)
        genome_list.append(genome)
        

        
    GOdistr=str(Counter(GOlist))[8:-1]
    GOset=','.join(str(s) for s in set(GOlist))        

    prod_distr=str(Counter(prod_list))[8:-1]
    distribution=len(set(genome_list))
    
    if distribution==113:
        OGclass="universal"
    elif distribution==1:
        OGclass='singleton'
    
    with_frameshift=frame_list.count("yes")
    all_hypothetical='no'
        
    if all([re.search("hypothetic",m) for m in prod_list]):
        product="hypothetical protein"
        all_hypothetical="yes"
    else:
        f=re.match("(.+?): [0-9]+",prod_distr[1:-1]).group(1)
        if re.search("hypothetical",f):
            product="hypothetical protein"
        else:
            product=f[1:-1]     

    out_list=[group, OGsize, distribution,OGclass, all_hypothetical,product,with_frameshift,prod_distr,trachomatis,psittaci,pneumoniae,pecorum,abortus,caviae,felis,muridarum]+counts+[GOdistr,GOset]
    out_list=[str(x) for x in out_list]
    out_list="\t".join(out_list)+"\n"
    out_f1.write(out_list) 
        
out_f1.close()

# add OG identifier to file with proteins info


try:
    os.remove(out1+'Chlamydia_proteins_v2.txt')
except:
    pass
out_f2 = open(out1+'Chlamydia_proteins_v2.txt', "a")

in_f5 = open(in1, "r")

tmp=in_f5.readline().rstrip()+"\tOGid\n"
out_f2.write(tmp)
for line in in_f5:
    line=line.rstrip()
    protein=line.split("\t")[0]
    out_line=line+"\t"+proteins_to_OG[protein]+"\n"
    out_f2.write(out_line)
    
out_f2.close()
    
