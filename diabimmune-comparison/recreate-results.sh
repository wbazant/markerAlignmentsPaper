# nextflow pull wbazant/marker-alignments-nextflow
# default params everywhere
# nextflow.config: 
# params {
#  ...
#  inputPath = 'DIABIMMUNE_WGS-sampleToFastqs.tsv'
#  libraryLayout = 'paired'
#  downloadMethod = 'wget'
#  unpackMethod = ''
#  bowtie2Command = 'bowtie2 --omit-sec-seq --no-discordant --no-unal -k10'
#  summarizeAlignmentsCommand = 'marker_alignments --min-read-query-length 60 --min-taxon-num-markers 2 --min-taxon-num-reads 2 --min-taxon-better-marker-cluster-averages-ratio 1.01 --threshold-avg-match-identity-to-call-known-taxon 0.97  --threshold-num-taxa-to-call-unknown-taxon 1 --threshold-num-markers-to-call-unknown-taxon 4     --threshold-num-reads-to-call-unknown-taxon 8'
# ...
# }


nextflow view wbazant/marker-alignments-nextflow | perl -nE 'say $1 if m{== content of file: (.*)}' | xargs -n1 dirname

for value in coverage cpm taxon_num_reads taxon_num_markers; do
  $(nextflow view wbazant/marker-alignments-nextflow | perl -nE 'say $1 if m{== content of file: (.*)}' | xargs -n1 dirname)/bin/makeTsv.pl results/summarizedAlignments .taxa.tsv $value \
    > diabimmune-our-results.$value.tsv
done

for value in coverage cpm taxon_num_reads taxon_num_markers; do 
  perl -E '
use strict;
use warnings;
use feature "say";
use List::MoreUtils qw/first_index/;
use File::Basename;


my ($dir, $pattern, $valueColumn) = @ARGV;
die "Usage: $0 dir pattern [column, default: "cpm"]" unless -d $dir && $pattern;
$valueColumn //= "cpm";

opendir(my $dh, $dir) || die "Can not opendir $dir: $!";
my @files = grep {$_=~m{$pattern}} readdir($dh);
closedir $dh;

die "No files in input directory" unless @files;

my %sampleLabels;

my %result;
for my $file (@files){
  open(my $fh, "<", "$dir/$file") or die "$!: $file";
  my $name = $file;
  if( -l "$dir/$file" ){
    my $f = readlink "$dir/$file";
    $name = basename $f;
  }
  my ($sample) = $name =~ m{(.*)$pattern};
  $sampleLabels{$sample}++;
  die $file unless $sample;
  my $valueColumnIndex;
  while(<$fh>){
    chomp;
    my @line = split "\t";

    if ($. == 1){
      $valueColumnIndex = first_index {$_ eq $valueColumn } @line;
      die "$file Header not recognised: $_" unless $valueColumnIndex > -1;
    } else{
      my $taxon = $line[0];
      my $value = $line[$valueColumnIndex];
      say join "\t", $sample, $taxon, $value;
    }
  }
}'  results/summarizedAlignments .taxa.tsv diabimmune-our-results-triples.$value.tsv
done

