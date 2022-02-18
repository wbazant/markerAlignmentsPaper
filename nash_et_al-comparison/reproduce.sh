# get "additional file 3" from https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0373-4#Sec20
if [ ! -f tmp/nash_et_al.xlsx ]; then
  curl "https://static-content.springer.com/esm/art%3A10.1186%2Fs40168-017-0373-4/MediaObjects/40168_2017_373_MOESM3_ESM.xlsx" > tmp/nash_et_al.xlsx
fi

#pip3 install openpyxl


if [ ! -f tmp/nash_et_al.tsv ]; then
  python -c 'import sys; import pandas; pandas.read_excel(sys.argv[1], sheet_name = "Table S3 - Fungal Reads", engine = "openpyxl").groupby(["SRA_Sample", "Name"]).size().to_csv(sys.stdout, sep = "\t")' tmp/nash_et_al.xlsx > tmp/nash_et_al.tsv
fi

if [ ! -f tmp/hmp.json ] ; then
  curl --globoff 'https://portal.hmpdacc.org/api/files?fields=file_format,file_type,file_annotation_pipeline,file_matrix_type&filters={}&from=1&save=&size=161265&sort=file_id:asc' > tmp/hmp.json
fi

if [ ! -f tmp/hmp.tsv ] ; then
  join -11 -21 -t $'\t' \
    <( cut -f1 tmp/nash_et_al.tsv | sort -u | tail -n+2  ) \
    <( jq -r '.data.hits | map(select(.file.subtype == "wgs_raw")) | map((.file.srs +"\t"+ .file.https))[]' < tmp/hmp.json | grep "v2" | sort -k1,1 -t $'\t') \
  > tmp/hmp.tsv
fi

python ./compare.py --nashXlsx "tmp/nash_et_al.xlsx" --ourTriplesTsv "taxon_num_reads.triples.tsv" --refdb-ncbi  "/home/wbazant/dev/markerAlignmentsPaper/refdb/taxa.sqlite" \
   > compare.out

