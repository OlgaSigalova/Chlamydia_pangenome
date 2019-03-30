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


query="""select p.*, f.nuc_seq_uncorrected from proteins p left join
        frameshifted_genes f
        on f.protein_id=p.protein_id
        where group_id='OG1000'"""
        
proteins=pd.read_sql(query, connection,index_col="protein_id")
proteins["seq"]=""

for i in proteins.index:
    if proteins.loc[i,"frameshift"]=="no":
        proteins.loc[i,"seq"]=proteins.loc[i,"nuc_seq"]
    else:
        proteins.loc[i,"seq"]=proteins.loc[i,"nuc_seq_uncorrected"]

df=pd.DataFrame(columns=["protein_id","N","from_start","to_end","length"])

n=8 # minimum length of polyN
m=0 # number of allowed mismatches in the stretch
for protein_id in proteins.index:
    p=proteins.loc[protein_id,"seq"]
    for N in ["A","C","G","T"]:
        (score,fine)=(0,0)
        from_end=0
        length=0
        l=len(p)
        for i in range(l):
            if p[i]==N:
                score+=1
            else:
                fine+=1
            if fine>m or i==l-1:
                if score>=n:
                    from_start=i-score-fine+1
                    to_end=l-i+score+fine-1
                    length=score
                    df=df.append(pd.DataFrame([[protein_id,N,from_start,to_end, length]],columns=["protein_id","N","from_start","to_end","length"]))
                (score,fine)=(0,0)


df.to_sql("polyN", con=connection, if_exists='replace',flavor='mysql',index=False)
