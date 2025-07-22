
### ==== 1. download data from TISCH database ====
python python tisch_data_download.py

### ==== 2. data selection based on QC and meta, prepare data for running, generate compass sbatch file ====
python tisch_data_prepare.py

### ==== 3. run COMPASS ====
for x in `ls tisch_mccp_compass`
do 
    sbatch -A cbp $x
done

### ==== 4. run Mebocost ===
h5_folder=./tisch_mccp_h5ad
for mccpid in `cut -f1 tisch2_selected_dataset_May28_2025_for_mccp.txt`
do
    ## === running
    echo $mccpid
    python tisch_mebo_run.py ${mccpid} ${h5_folder}/${mccpid}.h5ad
done
