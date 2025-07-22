import os,sys
import pandas as pd
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import h5py
import collections

init = 466 ## the last MCCP ID for existing datasets, the new data will be +1
tisch_celltype_dir='./tisch_celltype_ann/'
tisch_data_dir='./tisch_h5_data/'
h5ad_dir = './tisch_mccp_h5ad'
compass_dir = './tisch_mccp_compass/'


tisch_dat = pd.read_csv('./tisch2_selected_dataset_May28_2025.txt', sep = '\t', header = None)
tisch_dat.columns = ['TISCH_ID', 'dataset', 'species', 'therapy', 'patient_num', 'cell_num', 'platform', 'tumor', 'PMID']


## download the cell annotation file
for dataid in tisch_dat['dataset'].tolist():
    url = 'http://tisch.comp-genomics.org/static/data/%s/%s_CellMetainfo_table.tsv'%(dataid, dataid)
    os.system('wget -O %s/%s_CellMetainfo_table.tsv %s'%(tisch_celltype_dir, dataid, dataid))

## just realized some samples may have multiple conditions mixed within the sample
cond = ['Tissue', 'Site', 'Source', 'Treatment', 'Stage', 'Subtype', 'TNMstage']

cond_res = {}
for x in tisch_dat['dataset'].tolist():
    meta = pd.read_csv('data/%s_CellMetainfo_table.tsv'%x, sep = '\t')
    cond_res[x] = meta.loc[:,meta.columns.isin(cond)].apply(lambda col: collections.Counter(col.tolist()), axis = 0)


## currently to simplify, just keep those samples with one condition only
kept_dataset = []
for x in cond_res:
    c = []
    for y in cond_res[x].tolist():
        c.append(len(y)==1)
    if np.all(c) and len(cond_res[x]) == 1:
        kept_dataset.append(x)

kept_tisch_dat = tisch_dat[tisch_dat['dataset'].isin(kept_dataset)]

kept_tisch_dat['MCCP_ID'] = ['MCCP'+str(int(init)+x) for x in range(kept_tisch_dat.shape[0])]
kept_tisch_dat = kept_tisch_dat[['MCCP_ID']+[x for x in kept_tisch_dat.columns.tolist() if x != 'MCCP_ID']]

kept_tisch_dat.to_csv('./tisch2_selected_dataset_May28_2025_for_mccp_kept.txt', sep = '\t', index = None)


## download h5 file
for dataid in kept_tisch_dat['dataset'].tolist():
    url = 'http://tisch.comp-genomics.org/static/data/%s/%s_expression.h5'%(dataid, dataid)
    os.system('wget -O %s/%s_expression.h5 %s'%(tisch_data_dir, dataid, dataid))
