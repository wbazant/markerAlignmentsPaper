
ourResults() {
  join -11 -21 -t $'\t'  <( perl -nE 'my (@xs) = split "\t"; my ($gid) = m{20/(.*?)_R1}; say join "\t", $xs[0], $gid' DIABIMMUNE_WGS-sampleToFastqs.tsv  | sort -k1,1 ) <( cat ./diabimmune-our-results-triples.taxon_num_markers.tsv | sort -k1,1  )| perl -pe 's/\|/\t/' | cut -f 2,3 | sort | tr $'\t' : | tr -d ?
}
theirResults(){
  tail -n+2 diabimmune-results-published-with-eukdetect.csv | tr , $'\t' | cut -f 1,3 | sort | tr $'\t' :
}

allResultsTripleUsThemCommon(){
  comm <(ourResults ) <(theirResults)
}
commonResults(){
  comm -12 <(ourResults ) <(theirResults)
}

echo -n "Common results: "
commonResults | wc -l  

justUs(){
  comm -23 <(ourResults ) <(theirResults)
}

justThem(){
  comm -13 <(ourResults ) <(theirResults)
}


echo -n "Just us: "
justUs | wc -l

echo -n "Just them: "
justThem  | wc -l 
