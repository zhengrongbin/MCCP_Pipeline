import h5py
import numpy as np
import scipy.sparse as sp
import scanpy as sc
import pandas as pd

def _create_adata_(h5_path, meta_path):
    # 1) Open the file
    f = h5py.File(h5_path,'r')
    meta = pd.read_csv(meta_path, sep = '\t', index_col = 0)
    meta.index = meta.index.tolist()
    meta.columns = [x.replace('UMAP', 'umap') for x in meta.columns.tolist()]
    meta['cell_type'] = meta['Celltype (minor-lineage)'].tolist()
    
    # 2) Read the features table and pick out Gene Expression rows
    feat      = f['matrix/features']
    f_types   = feat['feature_type'][:]                 # array of bytes, e.g. b'Gene Expression'
    gene_idx  = np.where(f_types == b'Gene')[0]
    
    # 3) Grab the raw CSC arrays
    data    = f['matrix/data'][:]     # values
    indices = f['matrix/indices'][:]  # row indices for each value
    indptr  = f['matrix/indptr'][:]   # pointers into indices/data for each *column*
    shape   = tuple(f['matrix/shape'][:])  # (n_rows, n_cols)
    
    # 4) Build a CSC matrix (not CSR!)
    mat_csc = sp.csc_matrix((data, indices, indptr), shape=shape)
    
    # 5) Subset to Gene Expression *rows* and convert to CSR
    mat_csr = mat_csc.tocsr()[gene_idx, :]
    
    # 6) Decode barcodes & gene IDs
    barcodes = [bc.decode() for bc in f['matrix/barcodes'][:]]
    gene_ids = [g.decode()  for g in feat['id'][gene_idx]]
    
    # 7) Wrap in AnnData (cells × genes)
    adata = sc.AnnData(
        mat_csr.T,
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(index=gene_ids)
    )
    adata.obs = meta.reindex(index = adata.obs_names)
    ## sample down to 10,000 cells
    if adata.shape[0] > 10000:
        ## subsampling
        sc.pp.subsample(adata, n_obs = 10000, random_state = 123)
    
    return(adata)

tisch_celltype_dir='./tisch_celltype_ann/'
tisch_data_dir='./tisch_h5_data/'
h5ad_dir = './tisch_mccp_h5ad'
compass_dir = './tisch_mccp_compass/'

dat = pd.read_csv('tisch2_selected_dataset_May28_2025_for_mccp.txt', sep = '\t')
for i, line in dat.iterrows():
    Id = line['MCCP_ID']
    print(Id)
    dataset = line['dataset']
    exp_path = tisch_data_dir+'/%s_expression.h5'%(
        dataset
    )
    meta_path = tisch_celltype_dir+'/%s_CellMetainfo_table.tsv'%(
        dataset
    )
    adata = _create_adata_(exp_path, meta_path)
    adata.write(h5ad_dir+'/%s.h5ad'%Id)

### ===== average expression for compass =====
for i, line in dat.iterrows():
    Id = line['MCCP_ID']
    adata = sc.read_h5ad(h5ad_dir+'/%s.h5ad'%Id)
    n_celltype = adata.obs['cell_type'].unique().shape[0]
    ## set run time and RAM: 33 = 6hr, 1.53G; 8 = 1.5hr, 1.04G
    if n_celltype > 25:
        t = '16:00:00'
    elif n_celltype > 20:
        t = '12:00:00'
    elif n_celltype > 15:
        t = '8:00:00'
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
   
