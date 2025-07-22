import os,sys
import numpy as np
import pandas as pd
import collections
import pickle as pk
import cellxgene_census

init = 466 ## the last MCCP ID for existing datasets, the new data will be +1
cellxgene_data_dir='./cellxgene_h5ad_data/'
h5ad_dir = './cellxgene_mccp_h5ad'
compass_dir = './cellxgene_mccp_compass/'

## celltype meta
with cellxgene_census.open_soma(census_version="2025-01-30") as census:
    # Reads SOMADataFrame as a slice
    cell_metadata = census["census_data"]["homo_sapiens"].obs.read(
        column_names = ['dataset_id', "assay", "cell_type", "tissue", "tissue_general", "suspension_type", "disease"]
    )

    # Concatenates results to pyarrow.Table
    cell_metadata = cell_metadata.concat()

    # Converts to pandas.DataFrame
    cell_metadata = cell_metadata.to_pandas()


cell_metadata_tidy = cell_metadata.groupby(['dataset_id', 'tissue', 'disease']).apply(lambda df: df.iloc[0,].drop('cell_type').tolist()+[df.shape[0], df['cell_type'].unique().shape[0]])
cell_metadata_tidy = pd.DataFrame(cell_metadata_tidy.tolist(), columns = ['dataset_id', 'assay', 'tissue', 'tissue_general',
                                                    'suspension_type', 'disease', 'cell_num', 'celltype_num'])


## dataset meta
with cellxgene_census.open_soma(census_version="2025-01-30") as census:
    census_datasets = (
        census["census_info"]["datasets"]
        .read(column_names=["collection_name", "dataset_title", "dataset_id", "collection_doi", "collection_doi_label"])
        .concat()
        .to_pandas()
    )
    census_datasets = census_datasets.set_index("dataset_id")

## merge the two dataframe
dat_df = cell_metadata_tidy.merge(census_datasets, left_on = 'dataset_id', right_index = True)


## focus on popular platforms
dat_df = dat_df[dat_df['assay'].isin(["10x 3' v3", "10x 3' v2", "10x 5' v1", "Smart-seq2"])]

## ignore blood, cancer
dat_df = dat_df[~dat_df['tissue_general'].isin(['blood', 'lymph node']) &
                    ~dat_df['disease'].str.contains('carcinoma') &
                ~dat_df['disease'].isin(['metastatic melanoma', 'glioblastoma', 'B-cell non-Hodgkin lymphoma']) &
                ~dat_df['disease'].str.contains('tumor') &
                ~dat_df['disease'].str.contains('cancer')]
dat_df['tissue_disease'] = dat_df['tissue'].str.lower()+':'+dat_df['disease'].str.replace('normal', 'control').str.lower()

## exclude DISCO data to avoid redundant datasets
disco_df = pd.read_csv('disco_dataset_for_mccp_update.txt', sep = '\t')

disco_dataset = disco_df['tissue'].str.lower()+':'+disco_df['disease'].str.lower()
disco_dataset2 = disco_df['anatomical_site'].str.lower()+':'+disco_df['disease'].str.lower()

cg_dataset = dat_df['tissue_disease'].unique().tolist()

new = []
for c in cg_dataset:
    c = c.replace(' tissue', '')
    if c == 'skin of chest:control' or c == "skin of body:control":
        c = 'skin:normal'
    t, d = c.split(':')
    if d == "Crohn disease":
        d = "Crohn's disease"
    any_m = False
    mat_disco = [x for x in disco_dataset if d in x or x.split(':')[-1] in d]
    if not mat_disco:
        continue
    for dd in mat_disco:
        if t in dd or dd.split(':')[0] in t:
            any_m = True
    if not any_m:
        new.append(c)

new2 = []
for c in new:
    c = c.replace(' tissue', '')
    if c == 'skin of chest:control' or c == "skin of body:control":
        c = 'skin:normal'
    t, d = c.split(':')
    if d == "Crohn disease":
        d = "Crohn's disease"
    any_m = False
    mat_disco = [x for x in disco_dataset2 if d in x or x.split(':')[-1] in d]
    if not mat_disco:
        continue
    for dd in mat_disco:
        if t in dd or dd.split(':')[0] in t:
            any_m = True
    if not any_m:
        new2.append(c)

## new2 records new samples that are needed 
dat_df_select = dat_df[dat_df['tissue_disease'].str.lower().isin(new2)]

## cell QC filter
dat_df_select = dat_df_select.query('cell_num > 5000 and celltype_num > 5').sort_values('cell_num', ascending = False)
dat_df_select = dat_df_select[~dat_df_select['tissue_disease'].duplicated()]
dat_df_select = dat_df_select.sort_values('collection_name')
dat_df_select['MCCP_ID'] = ['MCCP'+str(int(init)+x) for x in range(dat_df_select.shape[0])]
dat_df_select = dat_df_select[['MCCP_ID']+[x for x in dat_df_select.columns.tolist() if x != 'MCCP_ID']]

dat_df_select.to_csv('cell_by_gene_selected_dataset_for_mccp_update.tsv', sep = '\t', index = None)

## download h5ad

Ids = dat_df_select['dataset_id'].unique().tolist()
for Id in Ids:
    print(Id)
    if os.path.exists(cellxgene_data_dir+'%s.h5ad'%Id):
        continue
    try:
        cellxgene_census.download_source_h5ad(dataset_id=Id, to_path=cellxgene_data_dir+'%s.h5ad'%Id)
    except:
        pass
        
