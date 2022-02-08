all: paper.pdf

refdb/ncbi_eukprot_met_arch_markers.fna:
	mkdir refdb
	cd refdb
	wget https://ndownloader.figshare.com/files/26173346
	tar -zxvf 26173346

refdbCrossValidation/nineTenth.fna.1.bt2: refdb/ncbi_eukprot_met_arch_markers.fna
	mkdir -pv refdbCrossValidation
	bash scripts/split_eukdb_reference.sh ./refdb/ ./refdbCrossValidation/oneTenth.fna ./refdbCrossValidation/oneTenthFolder ./refdbCrossValidation/nineTenth.fna
	bowtie2-build ./refdbCrossValidation/nineTenth.fna ./refdbCrossValidation/nineTenth.fna

refdbDoubled/doubled.fna.1.bt2: refdbCrossValidation/nineTenth.fna.1.bt2
	mkdir -pv refdbDoubled
	cat refdb/busco_taxid_link.txt > refdbDoubled/busco_taxid_link.txt
	cat refdb/busco_taxid_link.txt | perl -pe 's/at2759/bt2759/g'  >> refdbDoubled/busco_taxid_link.txt
	cat ./refdbCrossValidation/oneTenth.fna > refdbDoubled/doubled.fna
	cat ./refdbCrossValidation/oneTenth.fna | perl -pe 's/at2759/bt2759/g' | perl -pe 's/^/A/ unless /^>/'  >> refdbDoubled/doubled.fna
	bowtie2-build refdbDoubled/doubled.fna refdbDoubled/doubled.fna


unknownEuks/conf.yaml: refdbCrossValidation/nineTenth.fna.1.bt2
	mkdir -pv unknownEuks
	mkdir -pv unknownEuks/input
	mkdir -pv unknownEuks/results
	bash scripts/wgsim_one_tenth.sh refdbCrossValidation/oneTenthFolder refdbCrossValidation nineTenth.fna unknownEuks/input unknownEuks/conf.yaml unknownEuks/results

unknownEuks/results-summary.tsv: unknownEuks/conf.yaml
	perl -i -pe 's{bowtie2( --threads \d+)?( --seed \d+)?}{bowtie2 --threads 1 --seed 1337}g' ~/dev/EukDetect/rules/eukdetect.rules
	perl -i -pe 's{samtools view -q \d+ -bS}{samtools view -q 5 -bS}g' ~/dev/EukDetect/rules/eukdetect.rules
	snakemake --snakefile ~/dev/EukDetect/rules/eukdetect.rules --configfile unknownEuks/conf.yaml  --cores 1 runall
	head unknownEuks/results/*_hits_table.txt | perl -E '$$/ = "==>"; while(<>){my ($$header, @ls) = split "\n"; @ls = grep {not $$_ =~ m/^Name|^\w*$$|Empty read count|No taxa passing|==>/} @ls; my ($$fromSpecies) = $$header =~ m{results/(.*)_filtered_hits}; $$_ =~ s{\t.*}{} for @ls; say join "\t", $$fromSpecies, scalar @ls, @ls;  }' > unknownEuks/results-summary.tsv

unknownEuksBowtie2/results/our-method.results-summary.tsv: unknownEuks/conf.yaml
	mkdir -pv unknownEuksBowtie2
	mkdir -pv unknownEuksBowtie2/results.tmp
	mkdir -pv unknownEuksBowtie2/results
	ls `pwd`/unknownEuks/input/* | sort | perl -pe 's/\n/\t/ if $$. % 2 ' | perl -MFile::Basename -nE 'm{(.*).1.fq}; my $$x = basename $$1; print "$$x\t$$_"' > unknownEuksBowtie2/in.tsv
	bash scripts/run_our_method_on_unknown_euks.sh `pwd`/unknownEuksBowtie2/in.tsv `pwd`/refdbCrossValidation/nineTenth.fna `pwd`/refdb/busco_taxid_link.txt `pwd`/unknownEuksBowtie2/work `pwd`/unknownEuksBowtie2/results.tmp `pwd`/unknownEuksBowtie2/results

unknownEuksUnmodifiedEukdetect/results-summary.tsv: unknownEuks/results-summary.tsv
	mkdir -pv unknownEuksUnmodifiedEukdetect
	(cd unknownEuksUnmodifiedEukdetect && ln -sv ../unknownEuks/input input )
	mkdir -pv unknownEuksUnmodifiedEukdetect/results
	perl -pe 's/unknownEuks/unknownEuksUnmodifiedEukdetect/g' unknownEuks/conf.yaml > unknownEuksUnmodifiedEukdetect/conf.yaml
	perl -i -pe 's{bowtie2( --threads \d+)?( --seed \d+)?}{bowtie2 --threads 1 --seed 1337}g' ~/dev/EukDetect/rules/eukdetect.rules
	perl -i -pe 's{samtools view -q \d+ -bS}{samtools view -q 30 -bS}g' ~/dev/EukDetect/rules/eukdetect.rules
	
	snakemake --snakefile ~/dev/EukDetect/rules/eukdetect.rules --configfile unknownEuksUnmodifiedEukdetect/conf.yaml  --cores 1 runall
	head unknownEuksUnmodifiedEukdetect/results/*_hits_table.txt | perl -E '$$/ = "==>"; while(<>){my ($$header, @ls) = split "\n"; @ls = grep {not $$_ =~ m/^Name|^\w*$$|Empty read count|No taxa passing|==>/} @ls; my ($$fromSpecies) = $$header =~ m{results/(.*)_filtered_hits}; $$_ =~ s{\t.*}{} for @ls; say join "\t", $$fromSpecies, scalar @ls, @ls;  }' > unknownEuksUnmodifiedEukdetect/results-summary.tsv



supplement/wgsimWholeSamplesOneTenthCoverage.tsv: unknownEuks/results-summary.tsv unknownEuksBowtie2/results/our-method.results-summary.tsv unknownEuksUnmodifiedEukdetect/results-summary.tsv
	python3 scripts/parse_whole_samples_results.py \
		--refdb-marker-to-taxon-path refdb/busco_taxid_link.txt \
		--refdb-ncbi refdb/taxa.sqlite \
		--input "No filter:unknownEuksBowtie2/results/no-filter.results-summary.tsv" \
		--input "m2:unknownEuksBowtie2/results/m2.results-summary.tsv" \
		--input "M30:unknownEuksBowtie2/results/M30.results-summary.tsv" \
		--input "r4:unknownEuksBowtie2/results/r4.results-summary.tsv" \
		--input "m2M30:unknownEuksBowtie2/results/m2M30.results-summary.tsv" \
		--input "m2r4:unknownEuksBowtie2/results/m2r4.results-summary.tsv" \
		--input "r4M30:unknownEuksBowtie2/results/r4M30.results-summary.tsv" \
		--input "m2r4M30:unknownEuksBowtie2/results/m2r4M30.results-summary.tsv" \
		--input "EukDetect:unknownEuksUnmodifiedEukdetect/results-summary.tsv" \
		--input "modified EukDetect (MAPQ>=5):unknownEuks/results-summary.tsv" \
		--input "our method:unknownEuksBowtie2/results/our-method.results-summary.tsv" > supplement/wgsimWholeSamplesOneTenthCoverage.tsv

tmp:
	mkdir -pv tmp

tmp/wgsimMutationRate.json: tmp refdb/ncbi_eukprot_met_arch_markers.fna
	python3 scripts/simulate_and_align.py --refdb-ncbi refdb/taxa.sqlite --refdb-marker-to-taxon-path refdb/busco_taxid_link.txt --reference refdb/ncbi_eukprot_met_arch_markers.fna --sim-source refdb/ncbi_eukprot_met_arch_markers.fna --dir tmp --verbose --read-lengths 100 --base-error-rates 0.0 --mutation-rates 0.0 0.001 0.01 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2 > tmp/wgsimMutationRate.json

tmpDoubled:
	mkdir -pv tmpDoubled

tmpDoubled/wgsimMutationRate.json: tmpDoubled refdbDoubled/doubled.fna.1.bt2
	python3 scripts/simulate_and_align.py --refdb-ncbi refdb/taxa.sqlite --refdb-marker-to-taxon-path refdbDoubled/busco_taxid_link.txt --reference refdbDoubled/doubled.fna --sim-source refdbDoubled/doubled.fna --dir tmpDoubled --verbose --read-lengths 100 --base-error-rates 0.0 --mutation-rates 0.0 > tmpDoubled/wgsimMutationRate.json

tmpLeaveOneOut:
	mkdir -pv tmpLeaveOneOut

tmpLeaveOneOut/wgsimMutationRateLeaveOneOut.json: tmpLeaveOneOut refdbCrossValidation/nineTenth.fna.1.bt2
	python3 scripts/simulate_and_align.py --refdb-ncbi refdb/taxa.sqlite --refdb-marker-to-taxon-path refdb/busco_taxid_link.txt --reference refdbCrossValidation/nineTenth.fna --sim-source refdbCrossValidation/oneTenth.fna --dir tmpLeaveOneOut --verbose --read-lengths 100 --base-error-rates 0.0 --mutation-rates 0.0 0.01 0.025 0.05 0.075 0.1 0.125 > tmpLeaveOneOut/wgsimMutationRateLeaveOneOut.json

figures/wgsimMutationRate.png: tmp/wgsimMutationRate.json
	mkdir -pv figures
	python3 scripts/plot_recall_and_mapq_for_varying_mutation_rate.py --input-json tmp/wgsimMutationRate.json --output-png figures/wgsimMutationRate.png 

figures/valuesOverMutationRate.png: ./tmp/wgsimMutationRate.json
	python3 scripts/plot_over_mutation_rate.py --input-json ./tmp/wgsimMutationRate.json --output-png figures/valuesOverMutationRate.png --aggregation-level SameSpecies

figures/valuesOverMutationRateUnknownSpecies.png: tmpLeaveOneOut/wgsimMutationRateLeaveOneOut.json
	python3 scripts/plot_over_mutation_rate.py --input-json tmpLeaveOneOut/wgsimMutationRateLeaveOneOut.json --output-png figures/valuesOverMutationRateUnknownSpecies.png --aggregation-level SameGenus

figures/leaveOneOut.png: tmpLeaveOneOut/wgsimMutationRateLeaveOneOut.json
	mkdir -pv figures
	python3 scripts/plot_recall_and_mapq_for_varying_mutation_rate.py --input-json tmpLeaveOneOut/wgsimMutationRateLeaveOneOut.json --output-png figures/leaveOneOut.png

figures/bars.png: tmp/wgsimMutationRate.json
	mkdir -pv figures
	python3 scripts/plot_bars.py --input-json tmp/wgsimMutationRate.json --output-png figures/bars.png

figures/barsLeaveOneOut.png: tmpLeaveOneOut/wgsimMutationRateLeaveOneOut.json
	mkdir -pv figures
	python3 scripts/plot_bars.py --input-json tmpLeaveOneOut/wgsimMutationRateLeaveOneOut.json --output-png figures/barsLeaveOneOut.png 

figures/precisionBySpecies.png: tmp/wgsimMutationRate.json
	python3 scripts/plot_precision_by_species.py --input-alignments-sqlite tmp/100.0.0.0.0.alignments.sqlite --output-png figures/precisionBySpecies.png --refdb-ncbi refdb/taxa.sqlite --aggregation-level species

figures/dropoutForFilters.png: supplement/wgsimWholeSamplesOneTenthCoverage.tsv
	python3 scripts/plot_whole_samples_dropout_for_filters.py --input-tsv supplement/wgsimWholeSamplesOneTenthCoverage.tsv --output-png figures/dropoutForFilters.png

supplement/wgsim.tsv: tmp/wgsimMutationRate.json
	mkdir -pv supplement
	python3 scripts/wgsim_to_tsv.py  --input-json tmp/wgsimMutationRate.json --output-tsv supplement/wgsim.tsv

supplement/wgsimLeaveOneOut.tsv: tmpLeaveOneOut/wgsimMutationRateLeaveOneOut.json
	mkdir -pv supplement
	python3 scripts/wgsim_to_tsv.py  --input-json tmpLeaveOneOut/wgsimMutationRateLeaveOneOut.json --output-tsv supplement/wgsimLeaveOneOut.tsv

supplement/wgsimDoubled.tsv: tmpDoubled/wgsimMutationRate.json
	mkdir -pv supplement
	python3 scripts/wgsim_to_tsv.py  --input-json tmpDoubled/wgsimMutationRate.json --output-tsv supplement/wgsimDoubled.tsv

paper.pdf: paper.md biblio.bib figures/wgsimMutationRate.png figures/valuesOverMutationRate.png figures/valuesOverMutationRateUnknownSpecies.png  figures/leaveOneOut.png figures/bars.png figures/barsLeaveOneOut.png figures/precisionBySpecies.png supplement/wgsim.tsv  supplement/wgsimLeaveOneOut.tsv supplement/wgsimDoubled.tsv supplement/wgsimWholeSamplesOneTenthCoverage.tsv figures/dropoutForFilters.png
	pandoc -s --bibliography biblio.bib  --citeproc -f markdown paper.md -o paper.pdf
