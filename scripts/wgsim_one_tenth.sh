DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

oneTenthFolder="$1"
refdbDir="$2"
refdbPrefix="$3"
simPath="$4"
outputConf="$5"
outputResults="$6"

if ! ( [ "$oneTenthFolder" -a "$refdbDir" -a "$refdbPrefix" -a "$simPath" -a "$outputConf" -a "$outputResults" ] ) ; then
  echo "Usage: $0 oneTenthFolder refdbDir refdbPrefix simPath outputConf outputResults"
  exit 1
fi

#perl -i -pe 's{samtools view -q \d+ -bS}{samtools view -q 5 -bS}g' /home/wbazant/dev/EukDetect/rules/eukdetect.rules

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

for f in $oneTenthFolder/* ; do
  species=$(basename $f)
  echo "  $species:" >> $outputConf

  # num reads is number of letters divided by read length (100) multiplied by desired coverage (0.1)
  numReads=$(grep -v '>' ./refdbCrossValidation/oneTenthFolder/Yarrowia_keelungensis | perl -pe chomp | wc -c | perl -nE 'chomp; say sprintf("%d", $_ / 1000)' )
  wgsim -S 1337 -1100 -2100 -e 0.0 -r 0.0 -N $numReads $f $simPath/$species.1.fq  $simPath/$species.2.fq
done
