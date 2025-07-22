import os,sys
import collections
import pandas as pd
import numpy as np
from mebocost import mebocost
import h5py
import pickle as pk


meta_f = sys.argv[1]

combined_df_update = pd.read_csv(meta_f, sep = '\t')

## === 1. set comparison between disease vs normal with the same anatomic site, allowing different project 
can_comp = combined_df_update.groupby('Anatomic_Site').apply(lambda df: 'normal' in df['Disease'].tolist() and not np.all(df['Disease'] == 'normal'))
can_comp_datasets = can_comp[can_comp == True].index

comp_set1 = []
for x in can_comp_datasets.tolist():
    tmp = combined_df_update[combined_df_update['Anatomic_Site'] == x]
    normal = tmp.query('Disease == "normal"')
    for i, line in tmp.iterrows():
        if line['Disease'] != "normal":
            d = line['MCCP_ID']+'~'+line['Anatomic_Site']+'~'+line['Disease']+'~'+line['Projects']
            for n in (normal['MCCP_ID']+'~'+normal['Anatomic_Site']+'~'+normal['Disease']+'~'+normal['Projects']).tolist():
                comp_set1.append(d+'___vs___'+n)

### ==== 2. set comparison between disease and normal with the same project and tissue type
can_comp = combined_df_update.groupby(['Tissue', 'Projects']).apply(lambda df: 'normal' in df['Disease'].tolist() and not np.all(df['Disease'] == 'normal'))
can_comp_datasets = can_comp[can_comp == True].index

comp_set2 = []
for x, y in can_comp_datasets.tolist():
    tmp = combined_df_update[(combined_df_update['Tissue'] == x) & (combined_df_update['Projects'] == y)]
    normal = tmp.query('Disease == "normal"')
    for i, line in tmp.iterrows():
        if line['Disease'] != "normal":
            d = line['MCCP_ID']+'~'+line['Anatomic_Site']+'~'+line['Disease']+'~'+line['Projects']
            for n in (normal['MCCP_ID']+'~'+normal['Anatomic_Site']+'~'+normal['Disease']+'~'+normal['Projects']).tolist():
                comp_set2.append(d+'___vs___'+n)
            
## all comparison will be
comp_set = sorted(set(comp_set1+comp_set2))


### ==== 3. run diff mCCC using MEBOCOST =====
try: 
    os.mkdir('mebo_diff')
except:
    pass
    
for c in comp_set:
    cond1, cond2 = c.split('___vs___')
    mccp1 = cond1.split('~')[0]
    mccp2 = cond2.split('~')[0]
    ## load object
    obj1 = mebocost.load_obj(path = 'mebo_res/%s___mebov2_flux.pk'%mccp1)
    obj2 = mebocost.load_obj(path = 'mebo_res/%s___mebov2_flux.pk'%mccp2)
    ## concat
    obj = mebocost.concat_obj(obj1, obj2, cond1 = mccp1, cond2 = mccp2)
    ## run diff
    obj.CommDiff(comps=['%s_vs_%s'%(mccp1, mccp2)], thread = 8)
    if obj.diffcomm_res is None:
        os.system('echo +++++++%s'%c)
        continue
    if obj.diffcomm_res.get('%s_vs_%s'%(mccp1, mccp2), None) is None:
        os.system('echo +++++++%s'%c)
        continue
    obj.diffcomm_res['%s_vs_%s'%(mccp1, mccp2)].to_parquet('mebo_diff/%s_vs_%s___diffmCCC.parquet'%(mccp1, mccp2))
    os.system('echo +++++++%s'%c)

