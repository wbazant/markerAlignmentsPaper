DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
REF_PATH="~/dev/markerAlignmentsPaper/refdbCrossValidation"

nextflow run ./main.nf \
  --inputPath $DIR/in.tsv  \
  --resultDir $DIR/results \
  --downloadMethod local \
  --unpackMethod bz2 \
  --libraryLayout paired \
  --refdb ${REF_PATH}/nineTenth.fna \
  --markerToTaxonPath ${REF_PATH}/busco_taxid_link.txt  \
  --summarizeAlignmentsCommand "marker_alignments --min-taxon-num-markers 2 --min-taxon-num-reads 4" \
  -with-trace -resume | tee $DIR/tee.out
