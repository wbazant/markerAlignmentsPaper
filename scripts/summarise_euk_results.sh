#!/bin/bash

set -euo pipefail

if [ "$#" -eq 0 ]; then
   echo "Usage: $0 <dir>/*<filtered|all>_hits_table.txt"
   exit 1
fi

head -n1000 "$@" | perl -E '
$/ = "==>";
while(<>){
  my ($header, @ls) = split "\n";
  @ls = grep {not $_ =~ m/^Name|^\w*$|Empty read count|No taxa passing|==>/} @ls;
  my ($fromSpecies, @ms) = $header =~ m{.*/(.*)_(filtered|all)_hits};
  $_ =~ s{\t.*}{} for @ls;
  say join "\t", $fromSpecies, scalar @ls, @ls;
}
'
