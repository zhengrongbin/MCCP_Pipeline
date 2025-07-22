import os,sys
import time
from datetime import datetime
import scanpy as sc
import pandas as pd
import numpy as np
import argparse
from argparse import RawDescriptionHelpFormatter
from mebocost import mebocost
from scipy.sparse import csr_matrix
import scanpy as sc
import json
import pickle as pk

def compute_group_means(matrix, groups):
    # Convert all group labels to string to ensure homogeneity
    groups = np.array([str(g) for g in groups])
    unique_groups = np.unique(groups)
    group_means = {}
    
    for group in unique_groups:
        # Find indices for the current group
        group_indices = np.where(groups == group)[0]
        # Slice the matrix for these rows
        submatrix = matrix[group_indices, :]
        # Compute the mean along axis 0 (for each column)
        group_mean = np.squeeze(np.array(submatrix.mean(axis=0)))
        group_means[group] = group_mean
    
    return group_means


mccpid = sys.argv[1] ## unique dataset ID
h5_path = sys.argv[2] ## h5ad file

## conf file for mebocost
mebo_conf = '/lab-share/Cardio-Chen-e2/Public/rongbinzheng/MCCP/DISCO/NewDownload_March07_2025/mebocost.conf'
h5_folder = './disco_mccp_h5ad'
mebo_res_folder = './disco_mebo_res'
compass_dir = './disco_mccp_compass/'

if not os.path.exists(h5_path):
    sys.exit(0)

### ==== 1. read in h5ad ====
adata = sc.read_h5ad(h5_path)
adata = adata.raw.to_adata()
adata.obs_names = [adata.obs_names[i]+'_'+str(i) for i in range(adata.obs_names.shape[0])]

## === 2. run mebocost ====
mebo_obj = mebocost.create_obj(
                    adata = adata,
                    group_col = 'cell_type',
                    condition_col = None,
                    met_est = 'mebocost',
                    config_path = mebo_conf,
                    species='human',
                    cutoff_exp='auto',
                    cutoff_met='auto',
                    cutoff_prop=0.15,
                    sensor_type='All',
                    thread=8
                    )

## infer mCCC
mebo_obj.infer_commu(
                    n_shuffle=1000,
                    seed=123, 
                    Return=False, 
                    thread=None,
                    save_permuation=True,
                    min_cell_number = 10,
                    pval_method='permutation_test_fdr',
                    pval_cutoff=0.05
                    )

try:
    os.mkdir(mebo_res_folder)
except:
    pass
        
## save
mebocost.save_obj(obj=mebo_obj, path = mebo_res_folder+'/%s___mebov2.pk'%(mccpid))








