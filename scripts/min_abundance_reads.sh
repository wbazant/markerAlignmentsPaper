# adapted from wgsim_one_tenth
# takes first and last read
sourceFolder="$1"
refdbDir="$2"
refdbPrefix="$3"
simPath="$4"
outputConf="$5"
outputResults="$6"

if ! ( [ "$sourceFolder" -a "$refdbDir" -a "$refdbPrefix" -a "$simPath" -a "$outputConf" -a "$outputResults" ] ) ; then
  echo "Usage: $0 oneTenthFolder refdbDir refdbPrefix simPath outputConf outputResults"
  exit 1
fi

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
eukdetect_dir: "/home/wbazant/dev/EukDetect"
samples:
EOF

for f in $sourceFolder/*; do
  firstMarker=$( head -n4 $f | grep '^@' | cut -f 1 -d : | rev | cut -f 4- -d _ | rev )
  lastMarker=$( tail -n4 $f | grep '^@' | cut -f 1 -d : | rev | cut -f 4- -d _ | rev )
  if [ "$firstMarker" != "$lastMarker" ] ; then
    species=$(basename $f | cut -f 1 -d .)
    t=$simPath/$(basename $f)
    if [ $(basename $f) == "$species.1.fq" ] ; then
      echo "  $species:" >> $outputConf
    elif [ $(basename $f) != "$species.2.fq" ]; then
      echo "ERROR: $f"
      exit 1
    fi
    head -n4 $f > $t
    tail -n4 $f >> $t 
  fi    
done
