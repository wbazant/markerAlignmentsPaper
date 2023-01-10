set -euo pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
run_ma(){
  input="$1"
  nr="$2"
  output="$3"
  shift 1
  shift 1
  shift 1
  echo "$input -> $output"
  marker_alignments --refdb-marker-to-taxon-path $SCRIPT_DIR/../refdb/busco_taxid_link.txt --input $input  --output $output --output-type taxon_all --num-reads $nr $*
}
mkdir -pv $SCRIPT_DIR/{results-CORRAL,results-all,results-all-rl,results-all-rl-c}


for s in $SCRIPT_DIR/aln-all/*.a.sam; do 
  err=$(echo $s | rev | cut -f 1 -d / | rev | cut -f 1 -d .)
  num_reads=$(grep $err  $SCRIPT_DIR/read_counts.tsv | cut -f 2)

  output_corral=$SCRIPT_DIR/results-CORRAL/${err}.tsv
  test -f $output_corral || run_ma $s $num_reads $output_corral --min-read-query-length 60 --min-taxon-num-markers 2 --min-taxon-num-reads 2 --min-taxon-better-marker-cluster-averages-ratio 1.01 --threshold-avg-match-identity-to-call-known-taxon 0.97  --threshold-num-taxa-to-call-unknown-taxon 1 --threshold-num-markers-to-call-unknown-taxon 4     --threshold-num-reads-to-call-unknown-taxon 8 --num-reads $num_reads

  output_all=$SCRIPT_DIR/results-all/${err}.tsv
  test -f $output_all || run_ma $s $num_reads $output_all 

  output_rl=$SCRIPT_DIR/results-all-rl/${err}.tsv
  test -f $output_rl || run_ma $s $num_reads $output_rl --min-read-query-length 60

  output_rl_c=$SCRIPT_DIR/results-all-rl-c/${err}.tsv
  test -f $output_rl_c || run_ma $s $num_reads $output_rl_c --min-read-query-length 60 --min-taxon-num-markers 2 --min-taxon-num-reads 2 
done

mkdir -pv $SCRIPT_DIR/{results-best,results-best-rl,results-best-rl-c,results-best-rl-mapq,results-best-rl-mapq-c}
for s in $SCRIPT_DIR/aln-best/*.sam; do 
  err=$(echo $s | rev | cut -f 1 -d / | rev | cut -f 1 -d .)
  echo $err
  num_reads=$(grep $err  $SCRIPT_DIR/read_counts.tsv | cut -f 2)

  output_best=$SCRIPT_DIR/results-best/${err}.tsv
  test -f $output_best || run_ma $s $num_reads $output_best

  output_best_rl=$SCRIPT_DIR/results-best-rl/${err}.tsv
  test -f $output_best_rl || run_ma $s $num_reads $output_best_rl --min-read-query-length 60


  output_best_rl_c=$SCRIPT_DIR/results-best-rl-c/${err}.tsv
  test -f $output_best_rl_c || run_ma $s $num_reads $output_best_rl_c --min-read-query-length 60 --min-taxon-num-markers 2 --min-taxon-num-alignments 4

  output_best_rl_mapq=$SCRIPT_DIR/results-best-rl-mapq/${err}.tsv
  test -f $output_best_rl_mapq || run_ma $s $num_reads $output_best_rl_mapq --min-read-query-length 60 --min-read-mapq 30

  output_best_rl_mapq_c=$SCRIPT_DIR/results-best-rl-mapq-c/${err}.tsv
  test -f $output_best_rl_mapq_c || run_ma $s $num_reads $output_best_rl_mapq_c --min-read-query-length 60 --min-read-mapq 30 --min-taxon-num-markers 2 --min-taxon-num-alignments 4
done
