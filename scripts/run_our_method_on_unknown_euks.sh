#!/usr/bin/bash

inputPath="$1"
refdbPath="$2"
markerToTaxonPath="$3"
workDir="$4"
resultsTmp="$5"
resultsFinal="$6"

if ! ( [ "$1" -a "$2" -a "$3" -a "$4" -a "$5" -a "$6" ] ) ; then
  echo "Usage: $0 TODO"
  exit 1
fi

runOne(){
  outputPath="${resultsFinal}/$1"
  workDir="$2"
  bowtie2Command="$3"
  summarizeAlignmentsCommand="$4"
  if [ -f "$outputPath" ]; then
    return
  fi
  rm -fv "$resultsTmp/cpms.tsv"
  nextflow run wbazant/marker-alignments-nextflow -r main -w "$workDir" \
    --inputPath "$inputPath"  \
    --resultDir "$resultsTmp" \
    --downloadMethod local \
    --unpackMethod bz2 \
    --libraryLayout paired \
    --refdb "$refdbPath" \
    --markerToTaxonPath "$markerToTaxonPath" \
    --bowtie2Command "$bowtie2Command --threads 1" \
    --summarizeAlignmentsCommand "$summarizeAlignmentsCommand" \
    -with-trace -resume && python -c 'import sys; import pandas; df = pandas.read_csv(sys.stdin, sep="\t", header = 0, index_col = 0); [print(index + "\t" + str(col.notnull().sum()) + "\t" + "\t".join(col[col.notnull()].keys())) for index, col in df.iteritems()]' < "$resultsTmp/cpms.tsv" > "$outputPath"
}

runOne \
  no-filter.results-summary.tsv \
  "$workDir" \
  "bowtie2 --omit-sec-seq --no-discordant --no-unal" \
  "marker_alignments"

runOne \
  m2.results-summary.tsv \
  "$workDir" \
  "bowtie2 --omit-sec-seq --no-discordant --no-unal" \
  "marker_alignments --min-taxon-num-markers 2"

runOne \
  r4.results-summary.tsv \
  "$workDir" \
  "bowtie2 --omit-sec-seq --no-discordant --no-unal" \
  "marker_alignments  --min-taxon-num-reads 4"

runOne \
  M30.results-summary.tsv \
  "$workDir" \
  "bowtie2 --omit-sec-seq --no-discordant --no-unal" \
  "marker_alignments --min-read-mapq 30"

runOne \
  m2M30.results-summary.tsv \
  "$workDir" \
  "bowtie2 --omit-sec-seq --no-discordant --no-unal" \
  "marker_alignments --min-taxon-num-markers 2 --min-read-mapq 30"

runOne \
  r4M30.results-summary.tsv \
  "$workDir" \
  "bowtie2 --omit-sec-seq --no-discordant --no-unal" \
  "marker_alignments --min-taxon-num-reads 4 --min-read-mapq 30"

runOne \
  m2r4.results-summary.tsv \
  "$workDir" \
  "bowtie2 --omit-sec-seq --no-discordant --no-unal" \
  "marker_alignments --min-taxon-num-markers 2 --min-taxon-num-reads 4"

runOne \
  m2r4M30.results-summary.tsv \
  "$workDir" \
  "bowtie2 --omit-sec-seq --no-discordant --no-unal" \
  "marker_alignments --min-taxon-num-markers 2 --min-taxon-num-reads 4 --min-read-mapq 30"

runOne \
  our-method.results-summary.tsv \
  "$workDir" \
  "bowtie2 --omit-sec-seq --no-discordant --no-unal -k10" \
  "marker_alignments --min-read-query-length 60 --min-taxon-num-markers 2 --min-taxon-num-reads 2 --min-taxon-better-marker-cluster-averages-ratio 1.01 --threshold-avg-match-identity-to-call-known-taxon 0.97  --threshold-num-taxa-to-call-unknown-taxon 1 --threshold-num-markers-to-call-unknown-taxon 4     --threshold-num-reads-to-call-unknown-taxon 8"
