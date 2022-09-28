set -euo pipefail

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

inputConf="$1"
mapqFilter="$2"
eukdetectRulesFile="$3"
resultSummaryPath="$4"
resultSummaryUnfilteredPath="$5"

if ! ( [ "$inputConf" -a "$mapqFilter" -a "$eukdetectRulesFile" -a "$resultSummaryPath" -a "$resultSummaryUnfilteredPath" ] ) ; then
  echo "Usage: $0 inputConf mapqFilter eukdetectRulesFile resultSummaryPath resultSummaryUnfilteredPath"
  exit 1
fi

eukdetectResultDir=$(head $inputConf | perl -nE 'm{output_dir: "(.*)"} and print $1')
if [ ! "$eukdetectResultDir" ] ; then
  echo "$inputConf doesn't have output_dir at the top?"
  exit 1
fi

perl -i -pe 's{bowtie2( --threads \d+)?( --seed \d+)?}{bowtie2 --threads 1 --seed 1337}g' $eukdetectRulesFile
MAPQ_FILTER=$mapqFilter perl -i -pe 's{samtools view -q \d+ -bS}{samtools view -q $ENV{MAPQ_FILTER} -bS}g' $eukdetectRulesFile

mkdir -pv "$eukdetectResultDir"
mkdir -pv "$eukdetectResultDir/aln"
mkdir -pv "$eukdetectResultDir/filtering"
snakemake --snakefile $eukdetectRulesFile --configfile $inputConf  --cores 1 runall
bash $DIR/summarise_euk_results.sh $eukdetectResultDir/*filtered_hits_table.txt > $resultSummaryPath
bash $DIR/summarise_euk_results.sh $eukdetectResultDir/filtering/*all_hits_table.txt > $resultSummaryUnfilteredPath

