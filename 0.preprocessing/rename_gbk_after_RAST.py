import sys
import os
import fnmatch
import pandas as pd
import numpy as np

gbk_dir = "/home/bsgroup/Chlamydia/genomes_apr_2019"

# original names of gbk files
v1 = fnmatch.filter(os.listdir(gbk_dir), '*.gbk')
# new names (D + genome counter) 
v2 = ["D" + str(x + 1).zfill(3) + ".gbk" for x in range(len(v1))]


gbk_names_dict = dict(zip(v1, v2))

# save names dictionary
df = pd.DataFrame.from_dict(data = gbk_names_dict, orient='index')
df.columns = ["new_name"]
df.to_csv(gbk_dir + "/gbk_renamed_dict.csv")

# rename files
for f in v1:
    old_path = gbk_dir + "/" + f
    new_path = gbk_dir + "/" + gbk_names_dict[f]
    os.rename(old_path, new_path)
