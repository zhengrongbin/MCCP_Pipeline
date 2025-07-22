
### ==== 1. download data from cellxgene database ====
python cellxgene_data_download.py

### ==== 2. data selection based on QC and meta, prepare data for running, generate compass sbatch file ====
python cellxgene_data_prepare.py

### ==== 3. run COMPASS ====
for x in `ls cellxgene_mccp_compass`
do 
    sbatch -A cbp $x
done

### ==== 4. run Mebocost ===
h5_folder=./cellxgene_mccp_h5ad
for mccpid in `cut -f1 cell_by_gene_selected_dataset_for_mccp_update.tsv`
do
    ## === running
    echo $mccpid
    python cellxgene_mebo_run.py ${mccpid} ${h5_folder}/${mccpid}.h5ad
done
