import os,sys
import scanpy as sc
import pandas as pd
import numpy as np
import scanpy.external as sce

cellxgene_data_dir='./cellxgene_h5ad_data/'
h5ad_dir = './cellxgene_mccp_h5ad'
compass_dir = './cellxgene_mccp_compass/'



def _umap_(adata):
    adata.raw = adata.copy()
    adata = adata[:, adata.var.highly_variable]

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, random_state=123)

    sc.pp.neighbors(adata, random_state=123, n_pcs=30, n_neighbors=10)
    sc.tl.umap(adata, random_state=123)
    ## add umap cordinate to obs
    adata.obs['umap_1'] = adata.obsm['X_umap'][:,0]
    adata.obs['umap_2'] = adata.obsm['X_umap'][:,1]
    return(adata.raw.to_adata())
    

def _save_h5ad_(adata, mccp_id, to_dir):
    ## sample down to 10,000 cells
    if adata.shape[0] > 10000:
        ## subsampling
        sc.pp.subsample(adata, n_obs = 10000, random_state = 123)
    
    ## add umap cordinate to obs
    if 'X_umap' in adata.obsm:
        adata.obs['umap_1'] = adata.obsm['X_umap'][:,0]
        adata.obs['umap_2'] = adata.obsm['X_umap'][:,1]
    elif 'X_umap_scpoli' in adata.obsm:
        adata.obs['umap_1'] = adata.obsm['X_umap_scpoli'][:,0]
        adata.obs['umap_2'] = adata.obsm['X_umap_scpoli'][:,1]
    else:
        u = [x for x in adata.obsm if 'umap' in x]
        if len(u) > 0:
            adata.obs['umap_1'] = adata.obsm[u[0]][:,0]
            adata.obs['umap_2'] = adata.obsm[u[0]][:,1]
        else:
            adata = _umap_(adata)
    adata.write(os.path.join(to_dir, mccp_id+'.h5ad'))

def _create_h5ad_(dataset_df, h5_path, to_dir = './'):
    ## check dataset
    if dataset_df.shape[0] == 0:
        ## save to mccp id h5ad
        if os.path.exists(to_dir+'/'+dataset_df['MCCP_ID'].values[0]+'.h5ad'):
            return
        ## load h5ad
        adata = sc.read_h5ad(h5_path)
        ## del raw
        del adata.raw
        _save_h5ad_(adata=adata, mccp_id=dataset_df['MCCP_ID'].values[0], to_dir=to_dir)
    else:
        ## load h5ad
        adata = sc.read_h5ad(h5_path)
        ## del raw
        del adata.raw
        for i, line in dataset_df.iterrows():
            if os.path.exists(to_dir+'/'+line['MCCP_ID']+'.h5ad'):
                continue
            tissue, disease = line['tissue'], line['disease']
            tmp_adata = adata[(adata.obs['tissue'] == tissue) & (adata.obs['disease'] == disease)].copy()
            ## save to mccp id h5ad
            _save_h5ad_(adata=tmp_adata, mccp_id=line['MCCP_ID'], to_dir=to_dir)
    return

## ===== 1. prepare h5ad file for each dataset ====

samples = pd.read_csv('cell_by_gene_selected_dataset_for_updated_mccp.tsv', sep = '\t')
datasetid = samples['dataset_id'].unique()
for d in datasetid:
    tmp_df = samples.query('dataset_id == @d')
    h5_path = os.path.join(cellxgene_data_dir, d+'.h5ad')
    print(d)
    try:
        _create_h5ad_(tmp_df, h5_path, to_dir = h5ad_dir)
    except Exception as e:
        print(e)
        os.system('echo %s >> problem_sample.txt'%d)

## === 2. prepare average expression of cell type and provide sbatch for compass ===
mccp_all = os.listdir(h5ad_dir)

for s in mccp_all:
    Id = s.split('.')[0]
    adata = sc.read_h5ad('%s/%s'%(h5ad_dir,s))
    # adata.var_names = adata.var['feature_name'].copy()
    # adata.var_names_make_unique()
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
    avg_exp.to_csv('%s/%s.avg_exp.tsv'%(compass_dir, Id), sep = '\t')
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
    out = open('%s/%s_run.sbatch'%(compass_dir, Id), 'w')
    out.write(title)
    out.close()
    
