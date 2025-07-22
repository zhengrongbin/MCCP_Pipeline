
### ==== 1. download data from DISCO database ====
sh disco_data_download.sh

### ==== 2. data selection based on QC and meta, prepare data for running, generate compass sbatch file ====
python disco_data_collection.py

### ==== 3. run COMPASS ====
for x in `ls disco_mccp_compass`
do 
    sbatch -A cbp $x
done

### ==== 4. run Mebocost ===
h5_folder=./disco_mccp_h5ad
for mccpid in `cut -f1 disco_dataset_for_mccp_update.txt`
do
    ## === running
    echo $mccpid
    python disco_mebo_run.py ${mccpid} ${h5_folder}/${mccpid}.h5ad
done

