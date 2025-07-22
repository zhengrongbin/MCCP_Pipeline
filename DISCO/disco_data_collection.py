## 07-2025
## Author: Rongbin Zheng
## Note: the script can be run at E3 server in BCH, with required package installed,
## espically allowing running COMPASS using singularity, the compass_sh needs to be updated when any change applied

import os,sys
import pandas as pd
import numpy as np
import scanpy as sc
import scanpy.external as sce
from mebocost import mebocost
import collections
import matplotlib.pyplot as plt

init = 466 ## the last MCCP ID for existing datasets, the new data will be +1
disco_celltype_dir='./disco_celltype_ann/'
disco_data_dir='./disco_h5ad_data/'
h5ad_dir = './disco_mccp_h5ad'
compass_dir = './disco_mccp_compass/'
disco_meta = pd.read_csv('./metasample.txt', sep = '\t')

## ==== 1. data filter =====

## Focus on normal samples, and primary sample from disease biopsy
sample_type_need = ['control', 'disease tissue (non-cancer)'] 
## exclude samples, like cancer, cell line, because cancer cell annotation usually difficult, will use TISCH database
## cell line is not a good model for CCC analysis
disco_meta2 = disco_meta[disco_meta['sample_type'].isin(sample_type_need)]

## remove any genotypes, currently do not include samples with gene perturbations, such as KO, due to difficulties in interpration
disco_meta2 = disco_meta2[pd.isna(disco_meta2['genotype'])]

## remove sorting data, sorted cells will be difficult to model CCC due to one or two cell types available
disco_meta2 = disco_meta2[pd.isna(disco_meta2['cell_sorting'])]

## remove blood/fluid tissue, these cells usually flat and CCC analysis may be not that accurate as solid tissue
disco_meta2 = disco_meta2[~disco_meta2['tissue'].isin(['blood', 'umbilical cord blood', 'cerebrospinal fluid',
                                                     'blood/liver', 'urine', 'cord blood', 'pleural fluid',
                                                      'bone marrow and blood', 'peritoneal dialysis overnight effluent',
                                                      'bronchoalveolar lavage fluid'])]
disco_meta2['disease'] = [x if not pd.isna(x) else 'control' for x in disco_meta2['disease'].tolist()]

## remove treatment samples 
disco_meta2 = disco_meta2[pd.isna(disco_meta2['treatment']) | disco_meta2['treatment'].isin(['heart transplant'])]

### focus on high quality samples
## remove some samples with small cell numbers or umi numbers
disco_meta2 = disco_meta2.query('cell_number > 1000 and median_umi > 1000')
## collect all sample paths from folders
h5_pathes = {}
for i in range(1, 8):
    for x in os.listdir('%s/batch_%s/'%(disco_data_dir, i)):
        if x.endswith('.h5'):
            h5_pathes[x.replace('.h5', '')] = '%s/batch_%s/%s'%(disco_data_dir, i, x)


h5_path_col = []
problem_s = []
all_s = []
for i, line in df.iterrows():
    tmp = disco_meta2[(disco_meta2['project_id'] == line['project_id']) &
                (disco_meta2['tissue'] == line['tissue']) &
                (disco_meta2['disease'] == line['disease']) & 
                (disco_meta2['anatomical_site'] == line['anatomical_site'])]
    samples = tmp['sample_id'].tolist()
    all_s.extend(samples)
    ptmp = []
    for s in samples:
        if s in h5_pathes:
            ptmp.append(h5_pathes.get(s))
        else:
            problem_s.append(s)
    h5_path_col.append('; '.join(ptmp))

## sum cell number
df = disco_meta2.groupby(['project_id', 'tissue', 'disease', 'anatomical_site'])['cell_number'].sum().reset_index()

## set minimal cell number for a dataset
df_m = df[df['cell_number'] > 5000]
df_m = df_m.sort_values(['tissue', 'disease'])

## set up uniq ID
df_m['MCCP_ID'] = range(df_m.shape[0])
df_m['MCCP_ID'] = "MCCP"+df_m['MCCP_ID'].astype('str')


## cell type annotation 
## add cell type annotation file path
meta_path_col = []
for i, line in df_m.iterrows():
    tmp = disco_meta2[(disco_meta2['project_id'] == line['project_id']) &
                (disco_meta2['tissue'] == line['tissue']) &
                (disco_meta2['disease'] == line['disease']) & 
                (disco_meta2['anatomical_site'] == line['anatomical_site'])]
    ptmp = []
    for s in tmp['sample_id'].unique().tolist():
        path = '%s/%s.txt'%(disco_celltype_dir, s)
        if os.path.exists(path):
            ptmp.append(path)
        else:
            print(s)
    meta_path_col.append('; '.join(ptmp))
## attach celltype ann file path
df_m['meta_file'] = meta_path_col

### check cells with high confident cell type annotation cell_type_score > 0.6 as instructed by DISCO paper
cell_stat = {}
for i, line in df_m.iterrows():
    tmp_meta = pd.DataFrame()
    for f in line['meta_file'].split('; '):
        tmp_meta = pd.concat([tmp_meta, pd.read_csv(f, sep = '\t')])
    tmp_meta1 = tmp_meta.query('cell_type_score > 0.6')
    cell_stat[line['MCCP_ID']] = [tmp_meta.shape[0], tmp_meta['cell_type'].unique().shape[0], 
                                  tmp_meta1.shape[0], tmp_meta1['cell_type'].unique().shape[0],
                                  '; '.join(tmp_meta1['cell_type'].unique().tolist())]
    del tmp_meta
    del tmp_meta1

cell_stat = pd.DataFrame(cell_stat, index = ['orig_cell_ann_number', 'orig_cell_type_number',
                                'filter_cell_ann_number', 'filter_cell_type_number',
                                'filter_cell_type']).T
## merge cell stats into meta table
df_m = pd.merge(df_m, cell_stat, left_on = 'MCCP_ID', right_index = True)
## filter out low cell number and low cell type number
df_m = df_m.query('filter_cell_ann_number >= 5000 and filter_cell_type_number >= 5')

## re-assign unique ID
df_m = df_m.sort_values(['tissue', 'disease', 'anatomical_site'])
df_m['MCCP_ID'] = ['MCCP'+str(int(init)+x) for x in range(df_m.shape[0])]

## reorder
df_m = df_m[['MCCP_ID']+[x for x in df_m.columns.tolist() if x != 'MCCP_ID']]

## save out sample table
df_m.to_csv('disco_dataset_for_mccp_update.txt', sep = '\t', index = None)


## ====== 2. h5ad file and average expression preparation for COMPASS and mebocost =====

def _create_h5ad_(dataset_id, h5_path_list, meta_path_list, to_dir = './'):
    ## iterate meta and combine
    meta = pd.DataFrame()
    for m in meta_path_list:
        meta = pd.concat([meta, pd.read_csv(m, sep = '\t')])
    ## filter for cell_type_score
    meta = meta.query('cell_type_score > 0.6')

    ## iterate h5
    h5_list = []
    for h5 in h5_path_list:
        ## sample meta
        meta_tmp = meta[meta['sample_id'] == os.path.basename(h5).replace('.h5', '')]
        meta_tmp.index = meta_tmp['cell_id'].tolist()
        try:
            h5_tmp = sc.read_10x_h5(h5)
        except:
            os.system('echo %s >> problem_sample.txt'%h5)
            continue
        ## focus on meta cells
        h5_tmp = h5_tmp[h5_tmp.obs_names.isin(meta_tmp['cell_id'])]
        h5_tmp.obs = meta_tmp.loc[h5_tmp.obs_names]
        h5_list.append(h5_tmp)
        del h5_tmp
    adata = sc.concat(h5_list, axis = 0, join = 'outer')
    ## sample down to 10,000 cells
    if adata.shape[0] > 10000:
        ## subsampling
        sc.pp.subsample(adata, n_obs = 10000, random_state = 123)
    ## regular qc calculation
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    ## normalization
    sc.pp.normalize_total(adata, target_sum=1e04)
    # Logarithmize the data
    sc.pp.log1p(adata)
    ## try umap with fixed parameters, 
    ## the goal is to show cell types in umap plot, so it is ok without turning parameters for PCA and neighbors
    if adata.obs['sample_id'].unique().shape[0] > 1:
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample_id")
    else:
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)

    adata.raw = adata.copy()
    adata = adata[:, adata.var.highly_variable]

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, random_state=123)
    ## harmony to remove batch effects
    sce.pp.harmony_integrate(adata, "sample_id")
    adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']

    sc.pp.neighbors(adata, random_state=123, n_pcs=30, n_neighbors=10)
    sc.tl.umap(adata, random_state=123)
    ## add umap cordinate to obs
    adata.obs['umap_1'] = adata.obsm['X_umap'][:,0]
    adata.obs['umap_2'] = adata.obsm['X_umap'][:,1]
    adata.write(os.path.join(to_dir, dataset_id+'.h5ad'))

### ---- 2.1: prepare h5ad file -----
## sorted by total cell number and start with large one
samples = samples.sort_values('filter_cell_ann_number', ascending = False)

for i, line in samples.iterrows():
    dataset_id = line['MCCP_ID']
    if os.path.exists(h5ad_dir+dataset_id+'.h5ad'):
        continue
    h5_path_list = line['h5_path'].split('; ')
    meta_path_list = line['meta_file'].split('; ')
    print(dataset_id)
    try:
        _create_h5ad_(dataset_id, h5_path_list, meta_path_list, to_dir = h5ad_dir)
    except:
        os.system('echo %s >> problem_sample.txt'%dataset_id)

### ---- 2.2: average gene expression -----
mccp_all = [x for x in os.listdir(h5ad_dir) if x.endswith('.h5ad')]

for s in mccp_all:
    Id = s.split('.')[0]
    adata = sc.read_h5ad('%s/%s'%(h5ad_dir, s))
    adata = adata.raw.to_adata()
    n_celltype = adata.obs['cell_type'].unique().shape[0]
    ## set run time and RAM: 33 = 6hr, 1.53G; 8 = 1.5hr, 1.04G
    if n_celltype > 25:
        t = '20:00:00'
    elif n_celltype > 20:
        t = '15:00:00'
    elif n_celltype > 15:
        t = '10:00:00'
    else:
        t = '5:00:00'
    
    ### Running COMPASS for each cell type by the average gene expression
    ### output average gene expression
    avg_exp = sc.get.aggregate(adata, by = 'cell_type', func='mean')
    avg_exp = pd.DataFrame(avg_exp.layers['mean'], index = avg_exp.obs_names, columns = avg_exp.var_names).T
    # ## do un log since COMPASS will take log in the algorithm
    avg_exp = avg_exp.apply(lambda col: np.exp(col)-1)
    avg_exp.to_csv('%s/%s.avg_exp.tsv'%(compass_dir, Id, sep = '\t')
    ## generate sbatch for compass running
    title='''#!/bin/bash
#SBATCH --partition=bch-compute # queue to be used
#SBATCH --time=%s # Running time (in hours-minutes-seconds)
#SBATCH --job-name=%s # Job name
#SBATCH --mail-type=BEGIN,END,FAIL # send and email when the job begins, ends or fails
#SBATCH --mail-user=your_email_address # Email address to send the job status
#SBATCH --output=%s.txt # Name of the output file
#SBATCH --nodes=1 # Number of compute nodes
#SBATCH --ntasks=8 # Number of cpu cores on one node
#SBATCH --mem=3G

source /lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/bin/activate
module load singularity
echo "runing compass"
cond=%s
compass_dir=%s
exp_tsv=${compass_dir}/${cond}.avg_exp.tsv
out_dir=${compass_dir}/${cond}_compass_res
temp_dir=${compass_dir}/${cond}_temp

### this need to be updated if any changed for COMPASS software
singularity run --bind /lab-share/Cardio-Chen-e2/Public/rongbinzheng:/lab-share/Cardio-Chen-e2/Public/rongbinzheng /lab-share/Cardio-Chen-e2/Public/rongbinzheng/tmp/compass sh compass_run.sh ${exp_tsv} homo_sapiens ${out_dir} ${temp_dir} 8
'''%(t, Id, Id, Id, compass_dir)
    out = open('%s/%s_run.sbatch'%(Id, compass_dir), 'w')
    out.write(title)
    out.close()
    

