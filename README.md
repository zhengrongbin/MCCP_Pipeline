# A pipeline for MCCP data collection and processing 
MCCP stands for Metabolite-mediated Cell-cell Communication Portal. It integrates scRNA-seq atlas to map metabolite mediated cell-cell communications (mCCC) in different tissues, therefore, it facilitates new hypothesis generation and interesting discoveries for metabolite signaling using public scRNA-seq datasets. MCCP is available at http://rc-cbp.tch.harvard.edu/mccp.


## This pipline include stepwises instruction to collect and prepare dataset from the MCCP portal. The following sections were included:

- DISCO datasets: disco_pipeline.sh
- CellxGene datasets: cellxgene_pipeline.sh
- TISCH datasets: tisch_pipeline.sh
- Tidy-Up mCCC for individual dataset (e.g. integrate COMPASS result, generate UMAP plot and metabolite enzyme-sensor expression dot plot.)
- Perform differential mCCC between disease and normal datasets


## Installation for required packages and software
- Python packages: scanpy, pandas, numpy, matplotlib, seaborn, scipy, base64, h5py, cellxgene_census (https://github.com/chanzuckerberg/cellxgene-census), mebocost (https://github.com/kaifuchenlab/MEBOCOST). 
- Please install COMPASS refer to https://github.com/wagnerlab-berkeley/Compass.

## Detailed instruction for each database
### 1. DISCO datasets

- First, we use batch download API provided by DISCO database to download meta data and packed gene expression data in h5ad files. Please check the download links at https://disco.bii.a-star.edu.sg/download.
```{bash}
bash disco_data_download.sh
```
<p>After running the bash script, unzipped h5ad files will be stored into folder disco_h5ad_data. Cell annotation files will be stored into disco_celltype_ann. Meta information for all samples will be named as metasample.txt</p>

- Second, we select high quality samples with high quality cells based on the meta information and downloaded h5ad data. This selection was based on total cell number, cell type annotation accuracy, total UMI counts, etc.
```{bash}
python disco_data_collection.py

```
<p>By running this script, a new TXT file will be named as disco_dataset_for_mccp_update.txt for all samples that pass the filtering criteria, meaning high quality data. The script will also prepare h5ad file for running MEBOCOST and prepare averaged gene expression matrix of cell types for COMPASS running, as well as the average expression table and sbatch file to run compass (stored in folder of disco_mccp_compass).</p>

- Third, since COMPASS running take times, so can submit jobs to queue:
```{bash}
for x in `ls disco_mccp_compass`
do 
    sbatch -A cbp $x
done
```

- Fourth, we run MEBOCOST to infer mCCC for each scRNA-seq data.
```{bash}
h5_folder=./disco_mccp_h5ad
for mccpid in `cut -f1 disco_dataset_for_mccp_update.txt`
do
    ## === running
    echo $mccpid
    python disco_mebo_run.py ${mccpid} ${h5_folder}/${mccpid}.h5ad
done
```

- If the enough resource is available, the pipeline can be excute in a one-time manner:
```{bash}
bash disco_pipeline.sh
```

### 2. CellxGene datasets

- First, we download CellxGene dataset using cellxgene_census package provided by the database.
```{bash}
python cellxgene_data_download.py
```
<p>By running this, the sample meta file will be stored and named as cell_by_gene_selected_dataset_for_mccp_update.tsv. It also downloads h5ad file of gene expression into folder cellxgene_h5ad_data</p>

- Second, prepare h5ad file for MEBOCOST and average expression matrix for COMPASS.
```{bash}
python cellxgene_data_prepare.py
```
<p>This running gives h5ad ready to run MEBOCOST. Those h5ad files will be in folder named cellxgene_mccp_h5ad. Meanwhile, script and average expression matrix will be in cellxgene_mccp_compass folder.</p>

- Third, run COMPASS for flux estimation.
```{bash}
for x in `ls cellxgene_mccp_compass`
do 
    sbatch -A cbp $x
done
```

- Forth, run MEBOCOST for CellxGene datasets
```{bash}
h5_folder=./cellxgene_mccp_h5ad
for mccpid in `cut -f1 cell_by_gene_selected_dataset_for_mccp_update.tsv`
do
    ## === running
    echo $mccpid
    python cellxgene_mebo_run.py ${mccpid} ${h5_folder}/${mccpid}.h5ad
done
```


### 3. TISCH datasets

- First, copy the meta table of TISCH datasets at http://tisch.comp-genomics.org/gallery/, and save them into a TXT file named, such as, tisch2_selected_dataset_May28_2025.txt.
```{bash}
python tisch_data_download.py
```
<p>By running this, the sample meta file will be stored and named as tisch2_selected_dataset_May28_2025_for_mccp.txt. It also downloads h5 file of gene expression into folder tisch_h5_data and cell type annotation data into folder tisch_celltype_ann</p>

- Second, prepare h5ad file for MEBOCOST and average expression matrix for COMPASS.
```{bash}
python tisch_data_prepare.py
```
<p>This running gives h5ad ready to run MEBOCOST. Those h5ad files will be in folder named tisch_mccp_h5ad. Meanwhile, script and average expression matrix will be in tisch_mccp_compass folder.</p>

- Third, run COMPASS for flux estimation.
```{bash}
for x in `ls tisch_mccp_compass`
do 
    sbatch -A cbp $x
done
```

- Forth, run MEBOCOST for TISCH datasets
```{bash}
h5_folder=./tisch_mccp_h5ad
for mccpid in `cut -f1 tisch2_selected_dataset_May28_2025_for_mccp.txt`
do
    ## === running
    echo $mccpid
    python tisch_mebo_run.py ${mccpid} ${h5_folder}/${mccpid}.h5ad
done
```

### 4. Tidy up by integrating flux result and generating needed plots
<p>The following script will merge datasets from three repositories, and integate COMPASS flux to update mCCC results. The needed plots will be generated. All will be stored into mebo_res folder. It also produces a table containing all datasets, here, named as mccp_all_dataset_meta.tsv</p>

```{bash}
python integrate_compass_flux_and_prepare_data_plots.py disco_dataset_for_mccp_update.txt cell_by_gene_selected_dataset_for_mccp_update.tsv tisch2_selected_dataset_May28_2025_for_mccp.txt
```

### 5. Differential analysis between disease and normal

```{bash}
python differential_analysis.py mccp_all_dataset_meta.tsv
```
<p>Running the above script will generate comparisons between disease vs normal of samples from the same tissue types in the same project or the same anatomic site from different project. The differential mCCC result will be stored into folder named mebo_diff</p>








