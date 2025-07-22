### update mCCC tables using COMPASS result
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
## plot
import matplotlib
import seaborn as sns
from matplotlib import pyplot as plt
import pickle as pk
import io
import base64
import math
import scanpy as sc
import collections

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

def _umap_(df, title = ''):
    celltypes = df['cell_type'].unique().tolist()
    asg_id = pd.Series(range(0, len(celltypes))).astype('str')+': '+celltypes
    asg_id.index = celltypes
    df.insert(1, 'celltype', [asg_id[x] for x in df['cell_type'].tolist()])
    label_df = df.groupby(['celltype'])[['umap_1', 'umap_2']].mean()
    label_df.insert(1, 'label', [x.split(': ')[0] for x in label_df.index.tolist()])
    
    ncol = math.ceil(label_df.shape[0]/18)
    if ncol > 2:
        fig, ax = plt.subplots(figsize = (17, 5))
    else:
        fig, ax = plt.subplots(figsize = (12, 5))
    
    cm = sc.pl.palettes.godsnot_102
    NUM_COLORS = label_df.shape[0]
    
    node_col = {label_df.index.tolist()[i]:cm[i] for i in range(NUM_COLORS)}
    
    sns.scatterplot(data = df, x = 'umap_1', y = 'umap_2',
                   hue = 'celltype', edgecolor = 'none',
                   palette=node_col,#'tab20', 
                    s = 3, alpha = .6, zorder = 100)
    for i, line in label_df.iterrows():
        ax.text(line['umap_1'], line['umap_2'], line['label'], zorder = 100)
    sns.despine()
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), markerscale = 3, 
               ncol=math.ceil(label_df.shape[0]/18))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set(xlabel = 'UMAP_1', ylabel = 'UMAP_2', title = '')
    # plt.grid(zorder = -5)
    plt.tight_layout()
    ## save to bs4
    buffer4 = io.BytesIO()
    fig.savefig(buffer4, format='png', dpi = 180)
    buffer4.seek(0)
    png = base64.b64encode(buffer4.read()).decode('utf-8')
    plt.close()
    return(png)


def _dotplot_(plot_dat, x='index', c='avg_exp', s='perct',
           cmap = 'Purples', title = '', figsize = (7, 3.5)):
    fig, ax = plt.subplots(ncols=3, figsize=figsize, 
                           width_ratios = (5, .2, 1), 
                           gridspec_kw={'wspace': 0.05},
                          layout = 'constrained')
    
    sp = ax[0].scatter(
        x=plot_dat[x],
        y=[1] * plot_dat.shape[0],
        s=[0.5 if x == 0 else x for x in plot_dat[s]],
        c=plot_dat[c],
        cmap=cmap,
        linewidth = .5,
        edgecolors='black',
        vmax = np.percentile(plot_dat[c], 95)
    )
    ax[0].tick_params(axis='x', rotation=90)
    ax[0].set_yticks([])
    ax[0].set_ylabel('')  # no meaningful y‐axis
    ax[0].set_xlim(-0.5, len(plot_dat['index'])-0.5)
    ax[0].set_title(title, size = 8)
    # ── Add a colorbar for “avg_exp” ──
    
    cbar = fig.colorbar(sp, ax=ax[1])
    cbar.ax.tick_params(labelsize=8) 
    cbar.set_label('Mean', rotation=270, labelpad=10, fontsize = 8)
    
    # ── Build a custom legend for “perct” (dot sizes) ──
    # Choose a few representative size values (e.g. min, median, max)
    size_values = [
        plot_dat[s].min(),
        np.median(plot_dat['perct']),
        plot_dat[s].max()
    ]
    # Create “invisible” scatter handles for the legend
    handles = [
        ax[2].scatter([], [], 
                   s= .5 if sv == 0 else sv, 
                   color='black', 
                   edgecolors='none')
        for sv in size_values
    ]
    labels = [f"{int(sv)}" for sv in size_values]
    
    legend = ax[2].legend(
        handles,
        labels,
        fontsize = 8,
        title_fontsize = 8,
        # bbox_to_anchor = (1.13, 0),
        scatterpoints=1,
        frameon=False,
        alignment = 'left',
        title='% of above zero',
        loc='center right'
    )
    title = legend.get_title()
    title.set_ha('left')         # align horizontally to the left of its (x,y) coordinate
    title.set_va('center')       # align vertically to the center of its (x,y) coordinate  
    ax[2].add_artist(legend)
    for a in ax[1:]:
        a.spines['right'].set_visible(False)
        a.spines['top'].set_visible(False)
        a.spines['left'].set_visible(False)
        a.spines['bottom'].set_visible(False)
        a.set_xticks([])
        a.set_yticks([])
    ## save to bs4
    buffer4 = io.BytesIO()
    fig.savefig(buffer4, format='png', dpi = 180)
    buffer4.seek(0)
    png = base64.b64encode(buffer4.read()).decode('utf-8')
    plt.close()
    return(png)


# meta_f = sys.argv[1] ## disco_meta_summary table path
disco_file = sys.argv[1]
cellxgene_file = sys.argv[2]
tisch_file = sys.argv[3]

## == noted to check and update path if needed ====
disco_dir = './disco_mebo_res'
disco_compass = './disco_mccp_compass'

cellxgene_dir = './cellxgene_mebo_res'
cellxgene_compass = './cellxgene_mccp_compass'

tisch_dir = './tisch_mebo_res'
tisch_compass = './tisch_mccp_compass'

met_sen_path = '/lab-share/Cardio-Chen-e2/Public/rongbinzheng/mebocost_db/human_met_sensor_update_May21_2025.tsv'


disco_data = pd.read_csv(disco_file, sep = '\t')
cg_data = pd.read_csv(cellxgene_file, sep = '\t')
tisch_data = pd.read_csv(tisch_file, sep = '\t')
tisch_data = pd.read_csv('../TISCH/tisch2_selected_dataset_May28_2025_for_mccp_kept.txt', sep = '\t')
tisch_data['project_id'] = [x.split('_')[1] for x in tisch_data['dataset'].tolist()]
tisch_data['disease'] = [x.split('_')[0] for x in tisch_data['dataset'].tolist()]
cancer_disease_name = {'Glioma_GSE141383': 'Glioblastoma Multiforme', 
                       'NPC_GSE150430': 'Nasopharyngeal Carcinoma', 
                       'OSCC_GSE172577': 'Oral Squamous Cell Carcinoma',
                       'KIRC_GSE171306': 'Kidney Renal Clear Cell Carcinoma', 
                       'OS_GSE162454': 'Osteosarcoma', 
                       'OV_GSE147082': 'Metastatic Ovarian Serous Cystadenocarcinoma', 
                       'SS_GSE131309_Smartseq2': 'Synovial Sarcoma', 
                       'GIST_GSE162115':'Gastrointestinal Stromal Tumor', 
                       'THCA_GSE148673':'Thyroid Carcinoma'}

cancer_organ_name = {'Glioma_GSE141383': 'Brain', 
                       'NPC_GSE150430': 'Nasopharyngeal', 
                       'OSCC_GSE172577': 'Oral',
                       'KIRC_GSE171306': 'Kidney', 
                       'OS_GSE162454': 'Bone', 
                       'OV_GSE147082': 'omentum', 
                       'SS_GSE131309_Smartseq2': 'joint', 
                       'GIST_GSE162115':'Stomach', 
                       'THCA_GSE148673':'Thyroid'}

tisch_data['Organ'] = [cancer_organ_name.get(x) for x in tisch_data['dataset'].tolist()]
tisch_data['Disease'] = [cancer_disease_name.get(x) for x in tisch_data['dataset'].tolist()]
tisch_data['tissue'] = 'Tumor'

## uniform the table: MCCP_ID, organ, tissue, disease, project_id
cols = ['MCCP_ID','Tissue','Anatomic_Site','Disease','Projects', 'Cell_num']
## DISCO
df1 = disco_data[['MCCP_ID', 'tissue', 'anatomical_site', 'disease', 'project_id', 'cell_number']]
df1.columns = cols
df1['Disease'] = ['normal' if x == 'control' else x for x in df1['Disease'].tolist()]
df1['Source'] = 'DISCO'

## CellxGene
df2 = cg_data[['MCCP_ID', 'tissue_general', 'tissue', 'disease', 'collection_doi', 'cell_num']]
df2.columns = cols
df2['Source'] = 'CellxGene'

## TISCH
df3 = tisch_data[['MCCP_ID', 'Organ', 'tissue', 'Disease', 'project_id', 'cell_num']]
df3.columns = cols
df3['Source'] = 'TISCH'

### combine datasets
combined_df = pd.concat([df1, df2, df3])
combined_df['Tissue'] = combined_df['Tissue'].str.lower()
combined_df['Anatomic_Site'] = combined_df['Anatomic_Site'].str.lower()
combined_df['Disease'] = combined_df['Disease'].str.lower()
## some may don't have data
remove_datasets = []
for i, line in combined_df.iterrows():
    if line['Source'] == "DISCO":
        path = '%s/%s___mebov2.pk'%(disco_dir, line['MCCP_ID'])
    elif line['Source'] == "CellxGene":
        path = '%s/%s___mebov2.pk'%(cellxgene_dir, line['MCCP_ID'])
    elif line['Source'] == "TISCH":
        path = '%s/%s___mebov2.pk'%(tisch_compass, line['MCCP_ID'])
    else:
        print(i, 'Error')
        continue
    if not os.path.exists(path):
        remove_datasets.append(line['MCCP_ID'])

disco_data = combined_df[~combined_df['MCCP_ID'].isin(remove_datasets)]
disco_data.to_csv('mccp_all_dataset_meta.tsv', sep = '\t')


for i, line in disco_data.iterrows():
    if line['Source'] == "DISCO":
        path = os.path.join(disco_dir, '%s___mebov2.pk'%line['MCCP_ID'])
        compass_folder = os.path.join(disco_compass, '%s_compass_res/'%line['MCCP_ID'])
    elif line['Source'] == "CellxGene":
        path = os.path.join(cellxgene_dir, '%s___mebov2.pk'%line['MCCP_ID'])
        compass_folder = os.path.join(cellxgene_compass, '%s_compass_res/'%line['MCCP_ID'])
    elif line['Source'] == "TISCH":
        path = os.path.join(tisch_dir, '%s___mebov2.pk'%line['MCCP_ID'])
        compass_folder = os.path.join(tisch_compass, '%s_compass_res/'%line['MCCP_ID'])
    else:
        print(i, 'Error')
        continue
    
    mccpid = line['MCCP_ID']
    
    ## reload
    mebo_obj = mebocost.load_obj(path = path)
    ## update Sensor Annotation
    met_sen = pd.read_csv(met_sen_path, sep = '\t')
    sensor_type = {s: a for s, a in met_sen[['Gene_name', 'Annotation']].values.tolist()}
    
    ## for orignal table
    mebo_obj.original_result['Annotation'] = [sensor_type.get(x, np.nan) for x in mebo_obj.original_result['Sensor'].tolist()]
    ## comm table
    mebo_obj.commu_res['Annotation'] = [sensor_type.get(x, np.nan) for x in mebo_obj.commu_res['Sensor'].tolist()]
    
    ## apply constraint on compass flux result
    mebo_obj._ConstainCompassFlux_(compass_folder=compass_folder,
                                    efflux_cut='auto',
                                    influx_cut='auto',
                                    inplace=True)
    ### mCCC table to parquet
    mebo_obj.commu_res.to_csv('mebo_res/%s___mCCC.tsv.gz'%mccpid, sep = '\t', index = None, compression = 'gzip')
    mebocost.save_obj(mebo_obj, 'mebo_res/%s___mebov2_flux.pk'%mccpid)
    ### umap 
    df = mebo_obj.cell_ann[['umap_1', 'umap_2', 'cell_type']]
    umap_png = _umap_(df, title = mccpid)
    
    ## sensor expr 
    sensor_exp = collections.defaultdict()
    ## save the same celltype list
    uniq_celltype = mebo_obj.cell_ann['cell_type'].unique().tolist()
    sensor_exp['celltype'] = uniq_celltype
    for sensor in mebo_obj.commu_res['Sensor'].unique():
        ## get index of sensor expression from the sparse matrix
        indices = np.where(mebo_obj.exp_mat_indexer == sensor)[0]
        v = mebo_obj.exp_mat[indices,:]
        v = v.toarray()[0,]
        ## get celltype info as groups
        gr = mebo_obj.cell_ann.loc[mebo_obj.exp_mat_columns, 'cell_type']
        ## take avg expression based on group
        tmp_df = pd.DataFrame({'celltype': gr, 'expr': v})
        avg_exp = tmp_df.groupby('celltype')['expr'].mean()
        avg_exp = avg_exp.reindex(uniq_celltype).fillna(0)
        ## cal expression percentage (above 0) based on group
        perct = tmp_df.groupby('celltype').apply(lambda df: df.query('expr > 0').shape[0]*100/df.shape[0])
        perct = perct.reindex(uniq_celltype).fillna(0)
        ## save to a dict
        sensor_exp[sensor] = {'avg_exp': avg_exp.to_numpy(), 'perct': perct.to_numpy()}
    ## all sensor exp figures
    sensor_exp_fig = collections.defaultdict()
    for sensor in mebo_obj.commu_res['Sensor'].unique().tolist():
    
        plot_dat = pd.DataFrame(sensor_exp[sensor], index = sensor_exp['celltype']).reset_index()
        tl = np.max([len(x) for x in plot_dat['index'].tolist()])
        
        png = _dotplot_(plot_dat = plot_dat,
                  x='index', 
                  c='avg_exp', 
                  s='perct',
                  cmap = 'Reds', 
                  title = '%s expression' % sensor, 
                  figsize = (2.5+plot_dat.shape[0]*0.16, 1+tl*0.08))
        sensor_exp_fig[sensor] = png
        del png 
        del plot_dat
    
        
    ## met expr
    met_exp = collections.defaultdict()
    ## save the same celltype list
    uniq_celltype = mebo_obj.cell_ann['cell_type'].unique().tolist()
    met_exp['celltype'] = uniq_celltype
    for met in mebo_obj.commu_res['Metabolite'].unique():
        ## get index of sensor expression from the sparse matrix
        indices = np.where(mebo_obj.met_mat_indexer == met)[0]
        v = mebo_obj.met_mat[indices,:]
        v = v.toarray()[0,]
        ## get celltype info as groups
        gr = mebo_obj.cell_ann.loc[mebo_obj.met_mat_columns, 'cell_type']
        ## take avg expression based on group
        tmp_df = pd.DataFrame({'celltype': gr, 'expr': v})
        avg_exp = tmp_df.groupby('celltype')['expr'].mean()
        avg_exp = avg_exp.reindex(uniq_celltype).fillna(0)
        ## cal expression percentage (above 0) based on group
        perct = tmp_df.groupby('celltype').apply(lambda df: df.query('expr > 0').shape[0]*100/df.shape[0])
        perct = perct.reindex(uniq_celltype).fillna(0)
        ## save to a dict
        met_exp[met] = {'avg_exp': avg_exp.to_numpy(), 'perct': perct.to_numpy()}
    ## all met enzyme exp figures
    met_exp_fig = collections.defaultdict()
    for met, metname in mebo_obj.commu_res[['Metabolite', 'Metabolite_Name']].drop_duplicates().values.tolist():
    
        plot_dat = pd.DataFrame(met_exp[met], index = met_exp['celltype']).reset_index()
        tl = np.max([len(x) for x in plot_dat['index'].tolist()])
        
        png = _dotplot_(plot_dat = plot_dat,
                  x='index', 
                  c='avg_exp', 
                  s='perct',
                  cmap = 'Purples', 
                  title = '%s enzyme expression' % metname, 
                  figsize = (2.5+plot_dat.shape[0]*0.16, 1+tl*0.08))
        met_exp_fig[sensor] = png
        del png 
        del plot_dat

    try: 
        os.mkdir('mebo_res')
    except:
        pass
    with open('mebo_res/%s__met_sensor_umap_figs.pk'%mccpid, 'wb') as f:
        pk.dump({'met_png': met_exp_fig, 'sensor_png': sensor_exp_fig, 'umap_png': umap_png}, f)
    with open('mebo_res/%s__met_sensor_umap_data.pk'%mccpid, 'wb') as f:
        pk.dump({'met_dat': met_exp, 'sensor_dat': sensor_exp, 'umap_png': umap_png}, f)

        

    








