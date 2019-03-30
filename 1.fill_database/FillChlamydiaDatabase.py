#!/usr/bin/python

"""
Latest update on 16.05.2016
@author: olga
"""


from Bio import SeqIO
import sys
import os
import pandas as pd
import numpy as np
import MySQLdb
from collections import Counter 
import re
import itertools
import MySQLdb


###########################################################################
############ Connect and create databases #################################
###########################################################################


connection = MySQLdb.connect(
                host = 'centraldogma',
                user = 'bsgroup',
                passwd = 'Rerty4Z')  # create the connection
                
cursor = connection.cursor()
connection.autocommit(True)
cursor.execute("USE OGroups") # select the database

#create genome table

cursor.execute('DROP TABLE IF EXISTS genomes;')
cursor.execute("""CREATE TABLE genomes (
genome_id VARCHAR(10) NOT NULL,
PRIMARY KEY(genome_id),
full_name TEXT,
gi INT NOT NULL,
ref VARCHAR(20) NOT NULL, 
strain TEXT,
species TEXT,
DNAlength INT DEFAULT NULL,
serovar TEXT,
specificity TEXT,
cluster TEXT,
submitter TEXT,
year_sequenced TEXT,
collected_by TEXT,
year_collected TEXT,
anomality TEXT,
comment TEXT,
publication TEXT,
contigs INT);
""")

sql_gen=""" INSERT INTO genomes (genome_id, full_name, gi, ref, strain, species, DNAlength, 
        serovar, specificity, cluster, submitter, year_sequenced, collected_by, year_collected , 
        anomality, comment, publication, contigs) 
        VALUES (%s, %s, %s, %s,%s,%s, %s, %s, %s,%s,%s, %s, %s, %s,%s,%s,%s, %s);"""

cursor.execute('DROP TABLE IF EXISTS proteins;')
cursor.execute("""CREATE TABLE proteins (
protein_id VARCHAR(10) NOT NULL PRIMARY KEY,
genome_id VARCHAR(4) NOT NULL, 
product TEXT DEFAULT NULL,
start INT DEFAULT NULL,
end INT DEFAULT NULL,
strand VARCHAR(2) DEFAULT NULL,
length INT DEFAULT NULL,
frameshift TEXT NOT NULL,
nuc_seq TEXT DEFAULT NULL,
aa_seq TEXT NOT NULL,
group_id VARCHAR(10) NOT NULL,
group_id_withOutgroup VARCHAR(10) NOT NULL,
RAST2GO TEXT,
RAST2EC TEXT) 
""");

sql_prot="""INSERT INTO proteins (protein_id, genome_id, product,start, end, strand, length,
                                frameshift,nuc_seq,aa_seq, group_id, group_id_withOutgroup, RAST2GO, RAST2EC) 
                                VALUES (%s, %s, %s, %s,%s,%s, %s, %s, %s,%s,%s, %s, %s, %s);"""


cursor.execute('DROP TABLE IF EXISTS frameshifts;')
cursor.execute("""CREATE TABLE frameshifts (
frameshift_id VARCHAR(10) NOT NULL,
PRIMARY KEY(frameshift_id),
protein_id VARCHAR(10) NOT NULL,
fs_number INT,
to_end TEXT,
gap_size TEXT);
""")
sql_fsh = "INSERT INTO frameshifts (frameshift_id, protein_id, fs_number, to_end, gap_size) VALUES (%s,%s, %s, %s,%s);"


cursor.execute('DROP TABLE IF EXISTS frameshifted_genes;')
cursor.execute("""CREATE TABLE frameshifted_genes (
protein_id VARCHAR(10) NOT NULL,
PRIMARY KEY(protein_id),
num_frameshifts INT,
to_end TEXT,
gap_size TEXT,
nuc_seq_uncorrected TEXT);
""")
sql_fsg = "INSERT INTO frameshifted_genes (protein_id, num_frameshifts, to_end, gap_size,nuc_seq_uncorrected) VALUES (%s, %s, %s,%s,%s);"


cursor.execute('DROP TABLE IF EXISTS pfam;')
cursor.execute("""CREATE TABLE pfam (
pfam_id VARCHAR(20) NOT NULL PRIMARY KEY,
protein_id TEXT NOT NULL,
accession TEXT NOT NULL,
description TEXT DEFAULT NULL,
start INT DEFAULT NULL,
end INT DEFAULT NULL,
pfam2go TEXT DEFAULT NULL)""");

sql_pfam="""INSERT INTO pfam (pfam_id, protein_id, accession, description,start, end, pfam2go) 
                                VALUES (%s, %s, %s, %s,%s,%s, %s);"""


cursor.execute('DROP TABLE IF EXISTS ogroups;')
cursor.execute("""CREATE TABLE ogroups (
group_id VARCHAR(10) PRIMARY KEY,
group_size INT NOT NULL,
distribution INT NOT NULL,
length FLOAT NOT NULL,
all_hypothetical VARCHAR(3) NOT NULL,
product TEXT NOT NULL,
with_frameshift INT DEFAULT NULL,
product_distr TEXT,
trachomatis INT DEFAULT NULL,
psittaci INT DEFAULT NULL,
pneumoniae INT DEFAULT NULL,
pecorum INT DEFAULT NULL,
abortus INT DEFAULT NULL,
caviae INT DEFAULT NULL,
felis INT DEFAULT NULL,
muridarum INT DEFAULT NULL,
avium INT DEFAULT NULL,
gallinacea INT DEFAULT NULL,
ibidis INT DEFAULT NULL,
pfams TEXT DEFAULT NULL,
pfam_distr TEXT DEFAULT NULL,
GOdistr_pfam TEXT DEFAULT NULL,
GO_pfam TEXT DEFAULT NULL,
GOdistr_rast TEXT DEFAULT NULL,
GO_rast TEXT DEFAULT NULL,
ECdistr_rast TEXT DEFAULT NULL,
EC_rast TEXT DEFAULT NULL);""")

sql_og_withoutOutgroup = """INSERT INTO ogroups (group_id, group_size, distribution, length, 
                    all_hypothetical,product, with_frameshift,product_distr,trachomatis,
                    psittaci,pneumoniae,pecorum,abortus,caviae,felis,muridarum,avium,gallinacea, ibidis,
                    pfam_distr,pfams, GOdistr_pfam, GO_pfam, GOdistr_rast,GO_rast,
                    ECdistr_rast,EC_rast) VALUES (%s, %s,%s, %s,%s,%s, %s, %s,%s,
                    %s, %s, %s,%s,%s, %s, %s,%s,%s,%s, %s, %s,%s,%s, %s,%s,%s, %s);"""


cursor.execute('DROP TABLE IF EXISTS ogroups_nodraft;')
cursor.execute("""CREATE TABLE ogroups_nodraft (
group_id VARCHAR(10) PRIMARY KEY,
group_size INT NOT NULL,
distribution INT NOT NULL,
length FLOAT NOT NULL,
all_hypothetical VARCHAR(3) NOT NULL,
product TEXT NOT NULL,
with_frameshift INT DEFAULT NULL,
product_distr TEXT,
trachomatis INT DEFAULT NULL,
psittaci INT DEFAULT NULL,
pneumoniae INT DEFAULT NULL,
pecorum INT DEFAULT NULL,
abortus INT DEFAULT NULL,
caviae INT DEFAULT NULL,
felis INT DEFAULT NULL,
muridarum INT DEFAULT NULL,
avium INT DEFAULT NULL,
gallinacea INT DEFAULT NULL,
ibidis INT DEFAULT NULL,
pfams TEXT DEFAULT NULL,
pfam_distr TEXT DEFAULT NULL,
GOdistr_pfam TEXT DEFAULT NULL,
GO_pfam TEXT DEFAULT NULL,
GOdistr_rast TEXT DEFAULT NULL,
GO_rast TEXT DEFAULT NULL,
ECdistr_rast TEXT DEFAULT NULL,
EC_rast TEXT DEFAULT NULL);""")

sql_og_nodraft = """INSERT INTO ogroups_nodraft (group_id, group_size, distribution, length, 
                    all_hypothetical,product, with_frameshift,product_distr,trachomatis,
                    psittaci,pneumoniae,pecorum,abortus,caviae,felis,muridarum,avium,gallinacea, ibidis,
                    pfam_distr,pfams, GOdistr_pfam, GO_pfam, GOdistr_rast,GO_rast,
                    ECdistr_rast,EC_rast) VALUES (%s, %s,%s, %s,%s,%s, %s, %s,%s,
                    %s, %s, %s,%s,%s, %s, %s,%s,%s,%s, %s, %s,%s,%s, %s,%s,%s, %s);"""
   
                    

cursor.execute('DROP TABLE IF EXISTS ogroups_withOutgroup;')
cursor.execute("""CREATE TABLE ogroups_withOutgroup (
group_id VARCHAR(10) PRIMARY KEY,
group_size INT NOT NULL,
distribution INT NOT NULL,
length FLOAT NOT NULL,
all_hypothetical VARCHAR(3) NOT NULL,
product TEXT NOT NULL,
with_frameshift INT DEFAULT NULL,
product_distr TEXT,
trachomatis INT DEFAULT NULL,
psittaci INT DEFAULT NULL,
pneumoniae INT DEFAULT NULL,
pecorum INT DEFAULT NULL,
abortus INT DEFAULT NULL,
caviae INT DEFAULT NULL,
felis INT DEFAULT NULL,
muridarum INT DEFAULT NULL,
avium INT DEFAULT NULL,
gallinacea INT DEFAULT NULL,
ibidis INT DEFAULT NULL,
Waddlia_chondrophila INT DEFAULT NULL,
pfams TEXT DEFAULT NULL,
pfam_distr TEXT DEFAULT NULL,
GOdistr_pfam TEXT DEFAULT NULL,
GO_pfam TEXT DEFAULT NULL,
GOdistr_rast TEXT DEFAULT NULL,
GO_rast TEXT DEFAULT NULL,
ECdistr_rast TEXT DEFAULT NULL,
EC_rast TEXT DEFAULT NULL);""")

sql_og_withOutgroup = """INSERT INTO ogroups_withOutgroup (group_id, group_size, distribution, length, 
                    all_hypothetical,product, with_frameshift,product_distr,trachomatis,
                    psittaci,pneumoniae,pecorum,abortus,caviae,felis,muridarum,avium,gallinacea, ibidis,Waddlia_chondrophila,
                    pfam_distr,pfams, GOdistr_pfam, GO_pfam, GOdistr_rast,GO_rast,
                    ECdistr_rast,EC_rast) VALUES (%s, %s,%s, %s,%s,%s, %s, %s,%s,
                    %s, %s, %s,%s,%s, %s, %s,%s,%s,%s, %s, %s,%s,%s, %s,%s,%s,%s, %s);"""


########################################################################
############# Input - genomes info, ogs ################################
########################################################################

main_dir="/home/bsgroup/Chlamydia/"

# subset of relevant genomes
genomes=pd.read_csv(main_dir+"info/genomes_info_jun_2016.csv")
genomes=genomes[genomes["specificity"]!='recombinant'].replace(np.nan,'', regex=True)
filtered=['C091','C086','C112','C113','C017','C135','C131'] # these genomes are either duplicates (strain or whole genome) or of bad assembly
genomes=genomes[~genomes["genome_id"].isin(filtered)]
genomes.index=genomes["genome_id"]

genome_ids_withOutgroup=genomes["genome_id"].tolist() 
genome_ids_withoutOutgroup=[x for x in genome_ids_withOutgroup if x!="W001"] 
genome_ids_nodraft=genomes[(genomes["contigs"]==1) & (genomes["genome_id"]!="W001")]["genome_id"].tolist()    # without outgroup, no draft genomes

# table of Orthologous groups ( the code allows to make analysis on a subsample of OrthoMCL output. here - genome_ids)
OGs_withoutOutgroup= open(main_dir+'orthoMCL_jun_2016/results/Ogroups_without_outgroup_with_singletons.txt', "r")
OGs_withOutgroup=open(main_dir+'orthoMCL_jun_2016/results/Ogroups_with_outgroup_with_singletons.txt', "r")
prot2og_withoutOutgroup=dict()
og2prot_withoutOutgroup=dict()

prot2og_withOutgroup=dict()
og2prot_withOutgroup=dict()

# without outgroup, no draft genomes
prot2og_nodraft=dict()
og2prot_nodraft=dict()

# only genus Chlamydia (128 genomes)     
for line in OGs_withoutOutgroup:
    [group,proteins]=line.rstrip().split(":")
    proteins=proteins.strip().split(" ")
    proteins=[x.split("|")[1] for x in proteins if x.split("|")[0] in genome_ids_withoutOutgroup]
    og2prot_withoutOutgroup[group]=proteins
    for prot in proteins:
        prot2og_withoutOutgroup[prot]=group

# full sample (129 genomes)
for line in OGs_withOutgroup:
    [group,proteins]=line.rstrip().split(":")
    proteins=proteins.strip().split(" ")
    proteins=[x.split("|")[1] for x in proteins if x.split("|")[0] in genome_ids_withOutgroup]
    og2prot_withOutgroup[group]=proteins
    for prot in proteins:
        prot2og_withOutgroup[prot]=group       

# only full genomes of genus Chlamydia (116 genomes)
for line in OGs_withoutOutgroup:
    [group,proteins]=line.rstrip().split(":")
    proteins=proteins.strip().split(" ")
    proteins=[x.split("|")[1] for x in proteins if x.split("|")[0] in genome_ids_nodraft]
    og2prot_nodraft[group]=proteins
    for prot in proteins:
        prot2og_nodraft[prot]=group         

###########################################################################
################### Fill databases ########################################
###########################################################################
        
            
#################### 1. Genomes ###########################################
print ("parsing genomes information\n\n")

for i in range(genomes.shape[0]):
    cursor.execute(sql_gen,tuple(genomes.iloc[i,]))


#################### 2-4 proteins and frameshifts #########################
print ("parsing proteins information\n\n")

r=0 # counter for genomes
for genome_id in genome_ids_withOutgroup:
    input_handle  = open(main_dir+'genomes_jun_2016/'+genome_id+'.gbk', "r")
    gb_records = list(SeqIO.parse(input_handle, "genbank")) # relevant for gbk files with multiple contigs
    #gb_record = SeqIO.read(input_handle, "genbank")
    print str(r), genome_id
    r+=1
    k=0 # counter for CDSs per genome 
    n=0 # counter for frameshifts per genome
    for gb_record in gb_records:
        for seq_feature in gb_record.features :
            if seq_feature.type=="CDS" :
                k+=1
                OrthMCL_prot_code=genome_id+"_"+str(k)
                assert len(seq_feature.qualifiers['translation'])==1 # check that translation is given on one line
                start=seq_feature.location.start.position
                end=seq_feature.location.end.position
                strand=seq_feature.location.strand
                product=seq_feature.qualifiers['product'][0]
                try:
                    group_id=prot2og_withoutOutgroup[OrthMCL_prot_code] # will not work for W001
                except:
                    group_id=''
                group_id_withOutgroup=prot2og_withOutgroup[OrthMCL_prot_code]
                nucl_seq=seq_feature.extract(gb_record.seq) # sequences after correction for frameshifts
                aa_seq=seq_feature.qualifiers['translation'][0]
                lenght=len(aa_seq)
                frameshift="no"
				############ parse frameshifts ##########################
                if seq_feature.location_operator=='join':
                    frameshift="yes"
                    num_frameshifts=len(seq_feature.sub_features)-1
                    gap_size_arr=[] # distance between joined fragments. if negative, the fragment overlap
                    to_end_arr=[] # number of nucleotides from frameshift to gene end
                    if strand==-1:
                        nucl_seq_uc=gb_record.seq[start:end].reverse_complement() # to get uncorrected sequences
                        for i in range(num_frameshifts):
                            n+=1
                            frameshift_id=genome_id+"_fs"+str(n)
                            fs_number=num_frameshifts-i
                            gap_size=str(seq_feature.sub_features[i+1].location.start-seq_feature.sub_features[i].location.end)
                            gap_size_arr.append(gap_size)
                            to_end=str(seq_feature.sub_features[i].location.end-seq_feature.sub_features[0].location.start)
                            to_end_arr.append(to_end) 
                            cursor.execute(sql_fsh,(frameshift_id, OrthMCL_prot_code,fs_number,to_end,gap_size)) 
                    else:
                        nucl_seq__uc=gb_record.seq[start:end]
                        for i in range(num_frameshifts):
                            n+=1
                            frameshift_id=genome_id+"_fs"+str(n)
                            fs_number=i+1
                            gap_size=str(seq_feature.sub_features[i+1].location.start-seq_feature.sub_features[i].location.end)
                            gap_size_arr.append(gap_size)
                            to_end=str(seq_feature.sub_features[num_frameshifts].location.end-seq_feature.sub_features[i+1].location.start)
                            to_end_arr.append(to_end)
                            cursor.execute(sql_fsh,(frameshift_id, OrthMCL_prot_code,fs_number,to_end,gap_size)) 
                    gap_size_arr=";".join(gap_size_arr)
                    to_end_arr=";".join(to_end_arr)
                    cursor.execute(sql_fsg,(OrthMCL_prot_code,num_frameshifts,to_end_arr,gap_size_arr,nucl_seq__uc))
				#get EC numbers 
                EC=""
                try:
                    EC=seq_feature.qualifiers["EC_number"]
                    EC=";".join(EC)
                except:
                    pass
				#get GOs
                GOs=""
                if len(seq_feature.qualifiers['db_xref'])>1:
                    GOs=seq_feature.qualifiers['db_xref'][1:]
                    GOs=";".join(GOs)
				# write proteins info to database
                cursor.execute(sql_prot,(OrthMCL_prot_code, genome_id, product, start, end,strand,lenght,frameshift, nucl_seq, aa_seq, group_id,group_id_withOutgroup, GOs,EC))

######### Pfams ##################################################

print ("parsing pfams info\n\n")

pfam=pd.read_csv(main_dir+"info/InterProScan_results/Chlamydia_interproscan_w_outgroup_pfam.tsv",sep="\t", names=["protein_id", "Sequence_MD5_digest", "length", "analysis", "accession", "description", "start", "stop", "score", "status", "date", "ipr_annotation", "ipr_description", "GO"],index_col=False)
ids=set([x.split("|")[1] for x in pfam["protein_id"].tolist() if x.split("|")[0] in genome_ids_withOutgroup])
pfam_dict=dict((zip(ids,[1]*len(ids))))

for i in range (pfam.shape[0]):
    p=pfam.loc[i,"protein_id"].split("|")[1]
    try:
        pfam["id_short"]=p
        num=pfam_dict[p]
        pfam_dict[p]+=1
        pfam_id=p+"."+str(num)
        acc=pfam.loc[i,"accession"]
        desc=pfam.loc[i,"description"]
        start=pfam_descr=pfam.loc[i,"start"]
        end=pfam.loc[i,"stop"]
        pfam2go=pfam.loc[i,"GO"]
        if pd.isnull(pfam2go):
            pfam2go=""
        cursor.execute(sql_pfam,(pfam_id,p,acc,desc,start,end,pfam2go))
    except:
        pass

#################### Orthologous groups ################################

print ("creating summary for OGs")
prot=pd.read_sql('select * from proteins',connection,index_col='protein_id')

# without outgroup, 128 genomes

for group in og2prot_withoutOutgroup.keys():
    proteins=og2prot_withoutOutgroup[group]
    if len(proteins)>0: # excluded proteins (if some genomes are to be excluded after OGs are created)
        OGsize=len(proteins)
        prod_list=[]
        frame_list=[]
        genome_list=[]
        GO_rast=[]
        EC_rast=[]
        trachomatis, psittaci, pneumoniae, pecorum, abortus, caviae, felis, muridarum,avium,gallinacea, ibidis=[0,0,0,0,0,0,0,0,0,0,0]
        product=''
        distribution=0
        avg_length=prot[prot.index.isin(proteins)]["length"].mean()


        for p in proteins:

            #products
            prod=prot.loc[p,"product"]
            prod_list.append(prod)

            #frameshifts
            frame=prot.loc[p,"frameshift"]
            frame_list.append(frame)

            #species and genomes
            genome=prot.loc[p,"genome_id"]
            species=genomes.loc[genome,"species"]
            exec(species + "+= 1")
            genome_list.append(genome)


            #GO-terms
            if not pd.isnull(prot.loc[p,"RAST2GO"]):
                GOs=prot.loc[p,"RAST2GO"].split(";")
                for GO in GOs:
                    GO_rast.append(GO)


            #EC numbers
            if not pd.isnull(prot.loc[p,"RAST2EC"]):
                ECs=prot.loc[p,"RAST2EC"].split(";")
                for EC in ECs:
                    EC_rast.append(EC)
        
        #summary statistics
        GOdistr_rast=str(Counter(GO_rast))[8:-1]
        GO_rast=';'.join(str(s) for s in set(GO_rast)) 
        ECdistr_rast=str(Counter(EC_rast))[8:-1]
        EC_rast=';'.join(str(s) for s in set(EC_rast)) 
        prod_distr=str(Counter(prod_list))[8:-1]
        distribution=len(set(genome_list))
        with_frameshift=frame_list.count("yes")
        all_hypothetical='no'

        if all([re.search("hypothetic",m) for m in prod_list]):
            product="hypothetical protein"
            all_hypothetical="yes"
        else:
            f=re.match("(.+?): [0-9]+",prod_distr[1:-1]).group(1) # finds the most common product in the group
            if re.search("hypothetical",f):
                product="hypothetical protein"
            else:
                product=f[1:-1]    
        
        #pfam info
        if len(proteins)>1:
            query="select * from pfam where protein_id in "+str(tuple(proteins))
        else:
            query="select * from pfam where protein_id='"+str(proteins[0])+"'"
        df=pd.read_sql(query, connection)
        pfams=df["accession"].tolist()
        pfam_distr=str(Counter(pfams))[8:-1]
        pfams=';'.join(str(s) for s in set(pfams)) 
        #pfam to GO
        GO_pfam=df["pfam2go"].tolist()
        GO_pfam=[x.split("|") for x in GO_pfam]
        GO_pfam=list(itertools.chain.from_iterable(GO_pfam))
        GOdistr_pfam=str(Counter(GO_pfam))[8:-1]
        GO_pfam=';'.join(str(s) for s in set(GO_pfam)) 

        #add to database
        out_list=[group, OGsize, distribution,avg_length, all_hypothetical,product,with_frameshift,prod_distr,trachomatis, psittaci, pneumoniae, pecorum, abortus, caviae, felis, muridarum,avium,gallinacea, ibidis, pfam_distr,pfams,GOdistr_pfam,GO_pfam, GOdistr_rast,GO_rast, ECdistr_rast,EC_rast]
        cursor.execute(sql_og_withoutOutgroup,(out_list))
        #out_list=[str(x) for x in out_list]
        #out_list="\t".join(out_list)+"\n"


# with outgroup, 129 genomes

for group in og2prot_withOutgroup.keys():
    proteins=og2prot_withOutgroup[group]
    if len(proteins)>0: # excluded proteins (if some genomes are to be excluded after OGs are created)
        OGsize=len(proteins)
        prod_list=[]
        frame_list=[]
        genome_list=[]
        GO_rast=[]
        EC_rast=[]
        trachomatis, psittaci, pneumoniae, pecorum, abortus, caviae, felis, muridarum,avium,gallinacea, ibidis, Waddlia_chondrophila = [0,0,0,0,0,0,0,0,0,0,0,0]
        product=''
        distribution=0
        avg_length=prot[prot.index.isin(proteins)]["length"].mean()


        for p in proteins:

            #products
            prod=prot.loc[p,"product"]
            prod_list.append(prod)

            #frameshifts
            frame=prot.loc[p,"frameshift"]
            frame_list.append(frame)

            #species and genomes
            genome=prot.loc[p,"genome_id"]
            species=genomes.loc[genome,"species"]
            exec(species + "+= 1")
            genome_list.append(genome)


            #GO-terms
            if not pd.isnull(prot.loc[p,"RAST2GO"]):
                GOs=prot.loc[p,"RAST2GO"].split(";")
                for GO in GOs:
                    GO_rast.append(GO)


            #EC numbers
            if not pd.isnull(prot.loc[p,"RAST2EC"]):
                ECs=prot.loc[p,"RAST2EC"].split(";")
                for EC in ECs:
                    EC_rast.append(EC)
        
        #summary statistics
        GOdistr_rast=str(Counter(GO_rast))[8:-1]
        GO_rast=';'.join(str(s) for s in set(GO_rast)) 
        ECdistr_rast=str(Counter(EC_rast))[8:-1]
        EC_rast=';'.join(str(s) for s in set(EC_rast)) 
        prod_distr=str(Counter(prod_list))[8:-1]
        distribution=len(set(genome_list))
        with_frameshift=frame_list.count("yes")
        all_hypothetical='no'

        if all([re.search("hypothetic",m) for m in prod_list]):
            product="hypothetical protein"
            all_hypothetical="yes"
        else:
            f=re.match("(.+?): [0-9]+",prod_distr[1:-1]).group(1) # finds the most common product in the group
            if re.search("hypothetical",f):
                product="hypothetical protein"
            else:
                product=f[1:-1]    
        
        #pfam info
        if len(proteins)>1:
            query="select * from pfam where protein_id in "+str(tuple(proteins))
        else:
            query="select * from pfam where protein_id='"+str(proteins[0])+"'"
        df=pd.read_sql(query, connection)
        pfams=df["accession"].tolist()
        pfam_distr=str(Counter(pfams))[8:-1]
        pfams=';'.join(str(s) for s in set(pfams)) 
        #pfam to GO
        GO_pfam=df["pfam2go"].tolist()
        GO_pfam=[x.split("|") for x in GO_pfam]
        GO_pfam=list(itertools.chain.from_iterable(GO_pfam))
        GOdistr_pfam=str(Counter(GO_pfam))[8:-1]
        GO_pfam=';'.join(str(s) for s in set(GO_pfam)) 

        #add to database
        out_list=[group, OGsize, distribution,avg_length, all_hypothetical,product,with_frameshift,prod_distr,trachomatis, psittaci, pneumoniae, pecorum, abortus, caviae, felis, muridarum,avium,gallinacea, ibidis,Waddlia_chondrophila, pfam_distr,pfams,GOdistr_pfam,GO_pfam, GOdistr_rast,GO_rast, ECdistr_rast,EC_rast]
        cursor.execute(sql_og_withOutgroup,(out_list))
        #out_list=[str(x) for x in out_list]
        #out_list="\t".join(out_list)+"\n"


# without outgroup, no drafts



for group in og2prot_nodraft.keys():
    proteins=og2prot_nodraft[group]
    if len(proteins)>0: # excluded proteins (if some genomes are to be excluded after OGs are created)
        OGsize=len(proteins)
        prod_list=[]
        frame_list=[]
        genome_list=[]
        GO_rast=[]
        EC_rast=[]
        trachomatis, psittaci, pneumoniae, pecorum, abortus, caviae, felis, muridarum,avium,gallinacea, ibidis, Waddlia_chondrophila = [0,0,0,0,0,0,0,0,0,0,0,0]
        product=''
        distribution=0
        avg_length=prot[prot.index.isin(proteins)]["length"].mean()
        for p in proteins:
            #products
            prod=prot.loc[p,"product"]
            prod_list.append(prod)
            #frameshifts
            frame=prot.loc[p,"frameshift"]
            frame_list.append(frame)
            #species and genomes
            genome=prot.loc[p,"genome_id"]
            species=genomes.loc[genome,"species"]
            exec(species + "+= 1")
            genome_list.append(genome)
            #GO-terms
            if not pd.isnull(prot.loc[p,"RAST2GO"]):
                GOs=prot.loc[p,"RAST2GO"].split(";")
                for GO in GOs:
                    GO_rast.append(GO)
            #EC numbers
            if not pd.isnull(prot.loc[p,"RAST2EC"]):
                ECs=prot.loc[p,"RAST2EC"].split(";")
                for EC in ECs:
                    EC_rast.append(EC)
        #summary statistics
        GOdistr_rast=str(Counter(GO_rast))[8:-1]
        GO_rast=';'.join(str(s) for s in set(GO_rast)) 
        ECdistr_rast=str(Counter(EC_rast))[8:-1]
        EC_rast=';'.join(str(s) for s in set(EC_rast)) 
        prod_distr=str(Counter(prod_list))[8:-1]
        distribution=len(set(genome_list))
        with_frameshift=frame_list.count("yes")
        all_hypothetical='no'
        if all([re.search("hypothetic",m) for m in prod_list]):
            product="hypothetical protein"
            all_hypothetical="yes"
        else:
            f=re.match("(.+?): [0-9]+",prod_distr[1:-1]).group(1) # finds the most common product in the group
            if re.search("hypothetical",f):
                product="hypothetical protein"
            else:
                product=f[1:-1]    
        #pfam info
        if len(proteins)>1:
            query="select * from pfam where protein_id in "+str(tuple(proteins))
        else:
            query="select * from pfam where protein_id='"+str(proteins[0])+"'"
        df=pd.read_sql(query, connection)
        pfams=df["accession"].tolist()
        pfam_distr=str(Counter(pfams))[8:-1]
        pfams=';'.join(str(s) for s in set(pfams)) 
        #pfam to GO
        GO_pfam=df["pfam2go"].tolist()
        GO_pfam=[x.split("|") for x in GO_pfam]
        GO_pfam=list(itertools.chain.from_iterable(GO_pfam))
        GOdistr_pfam=str(Counter(GO_pfam))[8:-1]
        GO_pfam=';'.join(str(s) for s in set(GO_pfam)) 
        #add to database
        out_list=[group, OGsize, distribution,avg_length, all_hypothetical,product,with_frameshift,prod_distr,trachomatis, psittaci, pneumoniae, pecorum, abortus, caviae, felis, muridarum,avium,gallinacea, ibidis,pfam_distr,pfams,GOdistr_pfam,GO_pfam, GOdistr_rast,GO_rast, ECdistr_rast,EC_rast]
        cursor.execute(sql_og_nodraft,(out_list))

