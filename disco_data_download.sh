
## check the DISCO database for new data.tar.gz
## Date: March-07-2025
wget -O batch_1.tar.gz https://zenodo.org/records/14159931/files/batch_1.tar.gz?download=1
wget -O batch_2.tar.gz https://zenodo.org/records/14160154/files/batch_2.tar.gz?download=1
wget -O batch_3.tar.gz https://zenodo.org/records/14160213/files/batch_3.tar.gz?download=1
wget -O batch_4.tar.gz https://zenodo.org/records/14160221/files/batch_4.tar.gz?download=1
wget -O batch_5.tar.gz https://zenodo.org/records/14160748/files/batch_5.tar.gz?download=1
wget -O batch_6.tar.gz https://zenodo.org/records/14160802/files/batch_6.tar.gz?download=1
wget -O batch_7.tar.gz https://zenodo.org/records/14166702/files/batch_7.tar.gz?download=1

## unzip folder
outdir=./disco_h5ad_data
mkdir -p outdir

for i in 2 3 4 5 6 7
do
tar -xzvf batch_${i}.tar.gz -C ${outdir}
done

## ===== download meta sample =====
wget -O metasample.txt https://disco.bii.a-star.edu.sg/disco_v3_api/toolkit/getSampleMetadata


## download cell type table =====
while IFS=$'\t' read -r s p rest; do
    url = 'https://disco.bii.a-star.edu.sg/disco_v3_api/download/getMetadataTxt/%s/%s'%(p, s)
    out = open('disco_ann_download.sh', 'a')
    out.write('wget %s -O disco_celltype_ann/%s.txt\n'%(url, s))
    out.close()
done < metasample.txt

bash disco_ann_download.sh