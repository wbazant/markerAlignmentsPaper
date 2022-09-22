DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

refdb="$1"
oneTenth="$2"
oneTenthFolder="$3"
nineTenth="$4"
oneTenthCsv="$5"
nineTenthCsv="$6"

if ! ( [ "$refdb" -a "$oneTenth" -a "$oneTenthFolder" -a "$nineTenth" -a "$oneTenthCsv" -a "$nineTenthCsv" ] ) ; then
  echo "Usage: $0 refdb oneTenth oneTenthFolder nineTenth oneTenthCsv nineTenthCsv"
  exit 1
fi


cut -f3,4 -d, $refdb/marker_genes_per_species.csv | tail -n+2 | perl -nE 'print unless $. %10 == 0 ' > ${nineTenthCsv}

cut -f3,4 -d, $refdb/marker_genes_per_species.csv | tail -n+2 | perl -nE 'print if $. %10 == 0 ' > ${oneTenthCsv}


test -f ${nineTenth} || perl -E '
my %speciesToKeep;
open( my $fh, "<", $ARGV[0]) or die;
while(<$fh>){
  chomp;
  s{,.*}{};
  $speciesToKeep{$_}++;
}


$/ = ">";
open( my $fh, "<", $ARGV[1]) or die;
while(<$fh>){
  my @xs = split "-";
  next unless $speciesToKeep{$xs[1]};
  chomp;

  print ">$_";
}

' ${nineTenthCsv} $refdb/ncbi_eukprot_met_arch_markers.fna \
 > ${nineTenth}

test -f ${oneTenth} || perl -E '
my %speciesToKeep;
open( my $fh, "<", $ARGV[0]) or die;
while(<$fh>){
  chomp;
  s{,.*}{};
  $speciesToKeep{$_}++;
}


$/ = ">";
open( my $fh, "<", $ARGV[1]) or die;
while(<$fh>){
  my @xs = split "-";
  next unless $speciesToKeep{$xs[1]};
  chomp;

  print ">$_";
}
' ${oneTenthCsv} $refdb/ncbi_eukprot_met_arch_markers.fna \
 > ${oneTenth}


mkdir -pv ${oneTenthFolder}
perl -E '
my %speciesToKeep;
open( my $fh, "<", $ARGV[0]) or die;
while(<$fh>){
  chomp;
  s{,.*}{};
  $speciesToKeep{$_}++;
}


$/ = ">";
open( my $fh, "<", $ARGV[1]) or die;
my $outfh;
my $currentS;

while(<$fh>){
  my @xs = split "-";
  my $s = $xs[1];
  next unless $speciesToKeep{$s};
  if ($currentS ne $s){
    close $outfh if $outfh;
    $currentS = $s;
    open($outfh, ">", "$ARGV[2]/$s") or die $s;
  }

  chomp;
  print $outfh ">$_";
}
close $outfh if $outfh;
' ${oneTenthCsv} $refdb/ncbi_eukprot_met_arch_markers.fna ${oneTenthFolder}

