#!/bin/bash
set -euo pipefail
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

refdbDir="$1"
confusablePairsTsv="$2"
intermediateSimDir="$3"
simPath="$4"
outputConf="$5"
outputResults="$6"
eukdetectDir="$7"
refdbPrefix="$8"
coverage="$9"

if ! ( [ "$refdbDir" -a "$confusablePairsTsv" -a "$intermediateSimDir"  -a "$simPath" -a "$outputConf" -a "$outputResults" -a "$eukdetectDir" -a "$refdbPrefix" -a "$coverage" ] ) ; then
  echo "Usage: $0 refdbDir confusablePairsTsv intermediateSimDir simPath outputConf outputResults eukdetectDir refdbPrefix coverage"
  exit 1
fi


perl -nE 'if($. ==1) {next}; my ($ta, $na, $tb, $nb, $ipb, $ipf)  = split "\t"; next unless $ipb eq 'True' or $ipf eq 'True'; say "${ta}\t${tb}"' $confusablePairsTsv \
  > $intermediateSimDir/picked-pairs.tsv

join -t , -11 -24 \
  <( cat $intermediateSimDir/picked-pairs.tsv | xargs -n1 echo | sort -u ) \
  <( tail -n+2 $refdbDir/marker_genes_per_species.csv | sort -k4,4 -t , ) \
  | cut -f 1,4 -d, > $intermediateSimDir/species-picked-for-pairs.csv

mkdir -pv $intermediateSimDir/singles

perl -E '
my %speciesToKeep;
open( my $fh, "<", $ARGV[0]) or die;
while(<$fh>){
  chomp;
  my ($taxid, $name) = split ",";
  $speciesToKeep{$name} = $taxid;
}


$/ = ">";
open( my $fh, "<", $ARGV[1]) or die;
my $outfh;
my $currentS;

while(<$fh>){
  # handle species with dashes fungi-Fusarium_virguliforme_Mont-1-1314980at2759-S1
  my @xs = split "-";
  shift @xs;
  pop @xs;
  pop @xs;
  my $s = join "-", @xs;
  # refdb weird case where this does not match
  $s = "Diaporthe_longicolla_MSPL_11-6" if $s eq "Diaporthe_longicolla_MSPL_10-6";
  my $taxid = $speciesToKeep{$s};
  next unless $taxid;
  if ($currentS ne $s){
    close $outfh if $outfh;
    $currentS = $s;
    open($outfh, ">", "$ARGV[2]/${taxid}.fa") or die $s;
  }

  chomp;
  print $outfh ">$_";
}
close $outfh if $outfh;
' $intermediateSimDir/species-picked-for-pairs.csv $refdbDir/ncbi_eukprot_met_arch_markers.fna $intermediateSimDir/singles

cat <<EOF > $outputConf
output_dir: "$outputResults"
paired_end: true
fwd_suffix: ".1.fq"
rev_suffix: ".2.fq"
se_suffix: ".fastq.gz"
readlen: 100
fq_dir: "$simPath"
database_dir: "$refdbDir"
database_prefix: "$refdbPrefix"
eukdetect_dir: "$eukdetectDir"
samples:
EOF


cat $intermediateSimDir/picked-pairs.tsv \
  | while read taxonA taxonB; do 
  sample="taxonA${taxonA}taxonB${taxonB}"
  sourceA="$intermediateSimDir/singles/${taxonA}.fa"
  sourceB="$intermediateSimDir/singles/${taxonB}.fa"

  # num reads is number of letters divided by read length (100) multiplied by desired coverage (0.1)
  numReadsA=$(grep -v '>' $sourceA | perl -pe chomp | wc -c | CVG=$coverage perl -nE 'chomp; say sprintf("%d", ($_ / 100) * $ENV{CVG})' )
  numReadsB=$(grep -v '>' $sourceB | perl -pe chomp | wc -c | CVG=$coverage perl -nE 'chomp; say sprintf("%d", ($_ / 100) * $ENV{CVG})' )
  if [ $numReadsA -gt 0 -a $numReadsB -gt 0 ] ; then
    echo "  $sample:" >> $outputConf
    if [ -f "$simPath/$sample.1.fq" ] ; then
      echo "File exists: $simPath/$sample.1.fq, skipping wgsim for $sample"
    else
      wgsim -S 1337 -1100 -2100 -e 0.0 -r 0.0 -N $numReadsA $sourceA $intermediateSimDir/tmpA.1.fq $intermediateSimDir/tmpA.2.fq
      wgsim -S 1337 -1100 -2100 -e 0.0 -r 0.0 -N $numReadsB $sourceB $intermediateSimDir/tmpB.1.fq $intermediateSimDir/tmpB.2.fq
      cat $intermediateSimDir/tmpA.1.fq $intermediateSimDir/tmpB.1.fq > "$simPath/$sample.1.fq"
      cat $intermediateSimDir/tmpA.2.fq $intermediateSimDir/tmpB.2.fq > "$simPath/$sample.2.fq"
    fi
  fi
done
rm -f $intermediateSimDir/tmpA.1.fq $intermediateSimDir/tmpA.2.fq $intermediateSimDir/tmpB.1.fq $intermediateSimDir/tmpB.2.fq
