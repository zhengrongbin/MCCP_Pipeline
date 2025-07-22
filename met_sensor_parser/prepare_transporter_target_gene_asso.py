import pandas as pd
import numpy as np
import pyarrow.feather as feather
import re

### gene correlation from ARCHS4, downloaded at Nov 23 2021
### from page at https://maayanlab.cloud/archs4/download.html

mt_pair = pd.read_csv('/Users/rongbinzheng/Documents/github/Metabolic_Communication_V1/software/data/mascot_db/metabolite_sensors_annotation.tsv', sep = '\t') 
mt_pair['gene'] = mt_pair['Gene_name'].apply(lambda x: x.split('[')[0])
all_sensor_human = mt_pair['gene'].unique().tolist()

## human and mouse homology
homology = pd.read_csv('/Users/rongbinzheng/Documents/github/Metabolic_Communication_V1/software/data/mascot_db/human_mouse_homology_gene_pair.csv')
all_sensor_mouse = homology[homology['human_gene'].isin(all_sensor_human)]['mouse_gene'].unique().tolist()

### gene correlation from ARCHS4, downloaded at Nov 23 2021
path = '/Users/rongbinzheng/Documents/CommonData/ARCHS4/mouse/mouse_correlation_archs4.f'
read_df = feather.read_feather(path)
genes = [x[1:] if len(re.findall('^X[0-9]', x))>0 else x for x in read_df.columns]
read_df.index = genes
read_df.columns = genes
## all upper, not common mouse gene format, so check and reverse
mouse_gene = pd.read_csv('/Users/rongbinzheng/Documents/CommonData/mm10/gencode.vM23.annotation.protein_coding.csv')
mouse_gene['gene_upper'] = mouse_gene['gene_name'].str.upper()
mouse_gene = mouse_gene[~mouse_gene['gene_upper'].duplicated()]
mouse_gene.index = mouse_gene['gene_upper'].tolist()

## commom genes 
read_df = read_df.loc[read_df.index.isin(mouse_gene.index),
						read_df.columns.isin(mouse_gene.index)]

read_df.index = mouse_gene.loc[read_df.index.tolist(),'gene_name'].tolist()
read_df.columns = mouse_gene.loc[read_df.columns.tolist(),'gene_name'].tolist()

all_sensor_mouse = list(set(all_sensor_mouse) & set(read_df.columns.tolist()))

## extract only transporter x all gene 
need_corr = read_df[all_sensor_mouse]
## fisher z 
z_fun = lambda r: .5*(np.log(1+r) - np.log(1-r))
need_corr_fisherz = need_corr.apply(lambda col: z_fun(col))
need_corr_fisherz.to_csv('./mouse_corr_sensor_all_gene_pearson_r_fisherz.csv')



## ============ human =========

# mt_pair = pd.read_csv('./data/human/metabolite_transporter/manually_check_scFEA.txt') 
# all_transporter = mt_pair['transporter'].tolist()


### gene correlation from ARCHS4, downloaded at Nov 23 2021
path = '/Users/rongbinzheng/Documents/CommonData/ARCHS4/human/human_correlation_archs4.f'
read_df = feather.read_feather(path)
genes = [x[1:] if len(re.findall('^X[0-9]', x))>0 else x for x in read_df.columns]
read_df.index = genes
read_df.columns = genes

all_sensor_human = list(set(all_sensor_human) & set(read_df.columns.tolist()))
## extract only transporter x all gene 
need_corr = read_df[all_sensor_human]

z_fun = lambda r: .5*(np.log(1+r) - np.log(1-r))
need_corr_fisherz = need_corr.apply(lambda col: z_fun(col))
need_corr_fisherz.to_csv('./human_corr_sensor_all_gene_pearson_r_fisherz.csv')


