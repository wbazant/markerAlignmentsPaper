all: paper.pdf

refdb/ncbi_eukprot_met_arch_markers.fna:
	mkdir refdb
	cd refdb
	wget https://ndownloader.figshare.com/files/26173346
	tar -zxvf 26173346
	mv -v refdb/marker_genes_per_species.csv refdb/marker_genes_per_species.csv.orig
	perl -pe 's{Acartia_fossae,#N/A}{Acartia_fossae,1453877}; s{Calanus_sinicus,#N/A}{Calanus_sinicus,114070}' refdb/marker_genes_per_species.csv.orig > refdb/marker_genes_per_species.csv
	dos2unix refdb/marker_genes_per_species.csv

EukDetect-ad9edf11f5b458f11386b8a4b7f0e70f7bd69c30/rules/eukdetect.rules:
	wget https://github.com/allind/EukDetect/archive/ad9edf11f5b458f11386b8a4b7f0e70f7bd69c30.zip
	unzip ad9edf11f5b458f11386b8a4b7f0e70f7bd69c30.zip

refdbCrossValidation/nineTenth.fna.1.bt2: refdb/ncbi_eukprot_met_arch_markers.fna
	mkdir -pv refdbCrossValidation
	bash scripts/split_eukdb_reference.sh ./refdb/ ./refdbCrossValidation/oneTenth.fna ./refdbCrossValidation/oneTenthFolder ./refdbCrossValidation/nineTenth.fna ./refdbCrossValidation/oneTenth.csv ./refdbCrossValidation/nineTenth.csv
	bowtie2-build ./refdbCrossValidation/nineTenth.fna ./refdbCrossValidation/nineTenth.fna

refdbCrossValidation/goodMatches.tsv: refdbCrossValidation/nineTenth.fna.1.bt2
	python scripts/taxonomize_reference_split.py \
		--refdb-ncbi refdb/taxa.sqlite \
		--holdout-species-path ./refdbCrossValidation/oneTenth.csv \
		--remaining-species-path ./refdbCrossValidation/nineTenth.csv \
		--out-tsv-path refdbCrossValidation/goodMatches.tsv

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

unknownEuks/results-summary.tsv: unknownEuks/conf.yaml EukDetect-ad9edf11f5b458f11386b8a4b7f0e70f7bd69c30/rules/eukdetect.rules
	bash scripts/run_eukdetect_and_summarise.sh unknownEuks/conf.yaml 5 EukDetect-ad9edf11f5b458f11386b8a4b7f0e70f7bd69c30/rules/eukdetect.rules unknownEuks/results-summary.tsv unknownEuks/results-summary-unfiltered.tsv

lowAbundanceEuks/conf.yaml: unknownEuks/conf.yaml
	mkdir -pv lowAbundanceEuks
	mkdir -pv lowAbundanceEuks/input
	mkdir -pv lowAbundanceEuks/results
	bash scripts/min_abundance_reads.sh unknownEuks/input refdb ncbi_eukprot_met_arch_markers.fna lowAbundanceEuks/input lowAbundanceEuks/conf.yaml lowAbundanceEuks/results

lowAbundanceEuks/results-summary.tsv: lowAbundanceEuks/conf.yaml EukDetect-ad9edf11f5b458f11386b8a4b7f0e70f7bd69c30/rules/eukdetect.rules
	bash scripts/run_eukdetect_and_summarise.sh lowAbundanceEuks/conf.yaml 30 EukDetect-ad9edf11f5b458f11386b8a4b7f0e70f7bd69c30/rules/eukdetect.rules lowAbundanceEuks/results-summary.tsv lowAbundanceEuks/results-summary-unfiltered.tsv

lowAbundanceEuksBowtie2/results/our-method.results-summary.tsv: lowAbundanceEuks/conf.yaml
	mkdir -pv lowAbundanceEuksBowtie2
	mkdir -pv lowAbundanceEuksBowtie2/results.tmp
	mkdir -pv lowAbundanceEuksBowtie2/results
	ls `pwd`/lowAbundanceEuks/input/*fq | sort | perl -pe 's/\n/\t/ if $$. % 2 ' | perl -MFile::Basename -nE 'm{(.*).1.fq}; my $$x = basename $$1; print "$$x\t$$_"' > lowAbundanceEuksBowtie2/in.tsv
	bash scripts/run_CORRAL_with_variants.sh `pwd`/lowAbundanceEuksBowtie2/in.tsv `pwd`/refdbCrossValidation/nineTenth.fna `pwd`/refdb/busco_taxid_link.txt `pwd`/lowAbundanceEuksBowtie2/work `pwd`/lowAbundanceEuksBowtie2/results.tmp `pwd`/lowAbundanceEuksBowtie2/results

unknownEuksBowtie2/results/our-method.results-summary.tsv: unknownEuks/conf.yaml
	mkdir -pv unknownEuksBowtie2
	mkdir -pv unknownEuksBowtie2/results.tmp
	mkdir -pv unknownEuksBowtie2/results
	ls `pwd`/unknownEuks/input/*fq | sort | perl -pe 's/\n/\t/ if $$. % 2 ' | perl -MFile::Basename -nE 'm{(.*).1.fq}; my $$x = basename $$1; print "$$x\t$$_"' > unknownEuksBowtie2/in.tsv
	bash scripts/run_CORRAL_with_variants.sh `pwd`/unknownEuksBowtie2/in.tsv `pwd`/refdbCrossValidation/nineTenth.fna `pwd`/refdb/busco_taxid_link.txt `pwd`/unknownEuksBowtie2/work `pwd`/unknownEuksBowtie2/results.tmp `pwd`/unknownEuksBowtie2/results

unknownEuksBowtie2/results/our-method-unambiguous-only.results-summary.tsv: unknownEuksBowtie2/results/our-method.results-summary.tsv
	perl -nE 'chomp; my ($$id, $$num, @xs) = split ("\t", $$_, -1); @xs = grep {$$_} @xs; die "$$num bad: $$_" unless $$num == @xs; @ys = grep {$$_ !~/^\?/ } @xs ;say join "\t", $$id, scalar @ys, @ys' unknownEuksBowtie2/results/our-method.results-summary.tsv > unknownEuksBowtie2/results/our-method-unambiguous-only.results-summary.tsv


unknownEuksUnmodifiedEukdetect/results-summary.tsv: unknownEuks/results-summary.tsv EukDetect-ad9edf11f5b458f11386b8a4b7f0e70f7bd69c30/rules/eukdetect.rules
	mkdir -pv unknownEuksUnmodifiedEukdetect
	(cd unknownEuksUnmodifiedEukdetect && rm -rf input && ln -sv ../unknownEuks/input input)
	mkdir -pv unknownEuksUnmodifiedEukdetect/results
	perl -pe 's/unknownEuks/unknownEuksUnmodifiedEukdetect/g' unknownEuks/conf.yaml > unknownEuksUnmodifiedEukdetect/conf.yaml
	bash scripts/run_eukdetect_and_summarise.sh unknownEuksUnmodifiedEukdetect/conf.yaml 30 EukDetect-ad9edf11f5b458f11386b8a4b7f0e70f7bd69c30/rules/eukdetect.rules unknownEuksUnmodifiedEukdetect/results-summary.tsv unknownEuksUnmodifiedEukdetect/results-summary-unfiltered.tsv



unknownEuksBowtie2/results-summary-all.tsv: refdbCrossValidation/goodMatches.tsv unknownEuks/results-summary.tsv unknownEuksBowtie2/results/our-method.results-summary.tsv unknownEuksUnmodifiedEukdetect/results-summary.tsv unknownEuksBowtie2/results/our-method-unambiguous-only.results-summary.tsv
	python3 scripts/parse_whole_samples_results.py \
		--refdb-marker-to-taxon-path refdb/busco_taxid_link.txt \
		--refdb-ncbi refdb/taxa.sqlite \
		--good-matches-path refdbCrossValidation/goodMatches.tsv \
		--input "No filter:unknownEuksBowtie2/results/no-filter.results-summary.tsv" \
		--input "Two markers:unknownEuksBowtie2/results/m2.results-summary.tsv" \
		--input "MAPQ >=30:unknownEuksBowtie2/results/M30.results-summary.tsv" \
		--input "Four reads:unknownEuksBowtie2/results/r4.results-summary.tsv" \
		--input "Two markers, MAPQ>=30:unknownEuksBowtie2/results/m2M30.results-summary.tsv" \
		--input "Two markers, four reads:unknownEuksBowtie2/results/m2r4.results-summary.tsv" \
		--input "Four reads, MAPQ>=30:unknownEuksBowtie2/results/r4M30.results-summary.tsv" \
		--input "Two markers, four reads, MAPQ>=30:unknownEuksBowtie2/results/m2r4M30.results-summary.tsv" \
		--input "EukDetect (MAPQ>=30):unknownEuksUnmodifiedEukdetect/results-summary.tsv" \
		--input "EukDetect sensitive (MAPQ>=30):unknownEuksUnmodifiedEukdetect/results-summary-unfiltered.tsv" \
		--input "EukDetect (MAPQ>=5):unknownEuks/results-summary.tsv" \
		--input "EukDetect sensitive (MAPQ>=5):unknownEuks/results-summary-unfiltered.tsv" \
		--input "CORRAL (all hits):unknownEuksBowtie2/results/our-method.results-summary.tsv" \
		--input "CORRAL (unambiguous hits only):unknownEuksBowtie2/results/our-method-unambiguous-only.results-summary.tsv" \
		--output-tsv unknownEuksBowtie2/results-summary-all.tsv \
		--output-xlsx supplement/wgsimWholeSamplesOneTenthCoverage.xlsx


tmp:
	mkdir -pv tmp

tmp/wgsimMutationRate.json: tmp refdb/ncbi_eukprot_met_arch_markers.fna
	python3 scripts/simulate_and_align.py --refdb-ncbi refdb/taxa.sqlite --refdb-marker-to-taxon-path refdb/busco_taxid_link.txt --reference refdb/ncbi_eukprot_met_arch_markers.fna --sim-source refdb/ncbi_eukprot_met_arch_markers.fna --dir tmp --verbose --read-lengths 100 --base-error-rates 0.0 --mutation-rates 0.0 0.001 0.01 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2 > tmp/wgsimMutationRate.json

tmp/precisionBySpecies.tsv: tmp/wgsimMutationRate.json
	python3 scripts/sqlite_to_precision_by_species_tsv.py --input-alignments-sqlite tmp/100.0.0.0.0.alignments.sqlite --refdb-ncbi refdb/taxa.sqlite --aggregation-level species --output-tsv tmp/precisionBySpecies.tsv

tmpDoubled:
	mkdir -pv tmpDoubled

tmpDoubled/wgsimMutationRate.json: tmpDoubled refdbDoubled/doubled.fna.1.bt2
	python3 scripts/simulate_and_align.py --refdb-ncbi refdb/taxa.sqlite --refdb-marker-to-taxon-path refdbDoubled/busco_taxid_link.txt --reference refdbDoubled/doubled.fna --sim-source refdbDoubled/doubled.fna --dir tmpDoubled --verbose --read-lengths 100 --base-error-rates 0.0 --mutation-rates 0.0 > tmpDoubled/wgsimMutationRate.json

tmpLeaveOneOut:
	mkdir -pv tmpLeaveOneOut

tmpLeaveOneOut/wgsimMutationRateLeaveOneOut.json: tmpLeaveOneOut refdbCrossValidation/nineTenth.fna.1.bt2
	python3 scripts/simulate_and_align.py --refdb-ncbi refdb/taxa.sqlite --refdb-marker-to-taxon-path refdb/busco_taxid_link.txt --reference refdbCrossValidation/nineTenth.fna --sim-source refdbCrossValidation/oneTenth.fna --dir tmpLeaveOneOut --verbose --read-lengths 100 --base-error-rates 0.0 --mutation-rates 0.0 0.01 0.025 0.05 0.075 0.1 0.125 > tmpLeaveOneOut/wgsimMutationRateLeaveOneOut.json

tmpLeaveOneOut/precisionBySpeciesLeaveOneOut.tsv: tmpLeaveOneOut/wgsimMutationRateLeaveOneOut.json
	python3 scripts/sqlite_to_precision_by_species_tsv.py --input-alignments-sqlite tmpLeaveOneOut/100.0.0.0.0.alignments.sqlite --refdb-ncbi refdb/taxa.sqlite --aggregation-level genus --output-tsv tmpLeaveOneOut/precisionBySpeciesLeaveOneOut.tsv

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

# deprecated - Dan made a better figure
figures/precisionBySpecies.png: tmp/wgsimMutationRate.json
	python3 scripts/plot_precision_by_species.py --input-alignments-sqlite tmp/100.0.0.0.0.alignments.sqlite --output-png figures/precisionBySpecies.png --refdb-ncbi refdb/taxa.sqlite --aggregation-level species --output-tsv supplement/precisionBySpecies.tsv

figures/dropoutForFilters.png: unknownEuksBowtie2/results-summary-all.tsv
	echo "Needs fixing!"
	exit 1
	python3 scripts/plot_whole_samples_dropout_for_filters.py --input-tsv unknownEuksBowtie2/results-summary-all.tsv --output-png figures/dropoutForFilters.png

# was interesting once
supplement/crossvalidation_results_joined_with_num_reads.tsv: refdbCrossValidation/nineTenth.fna.1.bt2 unknownEuksBowtie2/results-summary-all.tsv
	join -11 -21 -t $$'\t' <( join -11 -21 -t $$'\t' <( cat refdbCrossValidation/busco_taxid_link.txt | perl -nE 'chomp; my ($$x, $$taxid) = split "\t"; my ($$s) = $$x =~ m{\w+-(.*)-\d+at2759.*}; next unless $$s; say join "\t", $$s, $$taxid' | sort -u ) <( grep -c '^>' refdbCrossValidation/oneTenthFolder/* | rev | cut -f 1 -d / | rev | tr : $$'\t' ) | cut -f 2,3 | perl -pE 'if($$.==1){say "taxid\tnum_reads"}'  | sort -r ) <(  perl -pE 'if($$.==1){s/^/taxid\t/}; s/\|/\t/;' unknownEuksBowtie2/results-summary-all.tsv | sort -r ) > supplement/crossvalidation_results_joined_with_num_reads.tsv


supplement/wgsim.tsv: tmp/wgsimMutationRate.json
	mkdir -pv supplement
	python3 scripts/wgsim_to_tsv.py  --input-json tmp/wgsimMutationRate.json --output-tsv supplement/wgsim.tsv

supplement/wgsimLeaveOneOut.tsv: tmpLeaveOneOut/wgsimMutationRateLeaveOneOut.json
	mkdir -pv supplement
	python3 scripts/wgsim_to_tsv.py  --input-json tmpLeaveOneOut/wgsimMutationRateLeaveOneOut.json --output-tsv supplement/wgsimLeaveOneOut.tsv

supplement/wgsimDoubled.tsv: tmpDoubled/wgsimMutationRate.json
	mkdir -pv supplement
	python3 scripts/wgsim_to_tsv.py  --input-json tmpDoubled/wgsimMutationRate.json --output-tsv supplement/wgsimDoubled.tsv

supplement/simulatedReads.xlsx: tmp/wgsimMutationRate.json tmp/precisionBySpecies.tsv tmpLeaveOneOut/wgsimMutationRateLeaveOneOut.json tmpLeaveOneOut/precisionBySpeciesLeaveOneOut.tsv
	python3 scripts/make_simulated_whole_reads_spreadsheet.py \
		--input-tsv-for-per-species-tab "Whole reference:tmp/precisionBySpecies.tsv" \
		--input-tsv-for-per-species-tab "Hold-out to remaining:tmpLeaveOneOut/precisionBySpeciesLeaveOneOut.tsv" \
  	--input-json-for-stats-tab "Whole reference:tmp/wgsimMutationRate.json" \
		--input-json-for-stats-tab "Hold-out to remaining:tmpLeaveOneOut/wgsimMutationRateLeaveOneOut.json" \
		--stats-tab-name "Stats" \
		--output-xlsx supplement/simulatedReads.xlsx

supplement/diabimmune.xlsx: diabimmune-comparison
	perl -nE 'my (@xs) = split "\t"; my ($$gid) = m{20/(.*?)_R1}; say join "\t", $$xs[0], $$gid'  diabimmune-comparison/DIABIMMUNE_WGS-sampleToFastqs.tsv > tmp.diabimmune-comparison-sample-to-run.tsv
	python3 scripts/make_diabimmune_spreadsheet.py \
		--refdb-ncbi refdb/taxa.sqlite \
		--input-tsv-triples-our-num-markers "diabimmune-comparison/diabimmune-our-results-triples.taxon_num_markers.tsv" \
		--input-tsv-triples-our-num-reads "diabimmune-comparison/diabimmune-our-results-triples.taxon_num_reads.tsv" \
		--input-csv-eukdetect "diabimmune-comparison/diabimmune-results-published-with-eukdetect.csv" \
		--input-tsv-sample-to-run "tmp.diabimmune-comparison-sample-to-run.tsv" \
		--output-xlsx diabimmune.tmp.xlsx && mv diabimmune.tmp.xlsx supplement/diabimmune.xlsx	&& rm tmp.diabimmune-comparison-sample-to-run.tsv

supplement/microbiomedb.xlsx: ./microbiomedb_results/BONUS.eukdetect.lineage_abundance.tsv
	python3 scripts/make_microbiomedb_stats_spreadsheet.py \
		--input-results-dir ./microbiomedb_results --input-name-postfix .eukdetect.lineage_abundance.tsv \
		--output-xlsx supplement/microbiomedb.xlsx


paper.pdf: paper.md biblio.bib figures/wgsimMutationRate.png figures/valuesOverMutationRate.png figures/valuesOverMutationRateUnknownSpecies.png  figures/leaveOneOut.png figures/bars.png figures/barsLeaveOneOut.png figures/precisionBySpecies.png supplement/wgsim.tsv  supplement/wgsimLeaveOneOut.tsv supplement/wgsimDoubled.tsv unknownEuksBowtie2/results-summary-all.tsv supplement/simulatedReads.xlsx supplement/diabimmune.xlsx supplement/microbiomedb.xlsx
	perl -pe 's/≥/\$$\\geq\$$/g; s/μ/\$$\\mu\$$/g; s/≤/\$$\\leq\$$/g' paper.md >  out.md
	pandoc -s --bibliography biblio.bib  --citeproc --csl bmc-bioinformatics.csl -f markdown out.md  --pdf-engine=xelatex -o paper.pdf
	perl -pe 's/≥/>=/g; s/μ/M/g; s/≤/<=/g; s/for .*?lust.*?ignments/for Clustering Of Related Reference ALignments/g' paper.md >  out.md
	pandoc -s --bibliography biblio.bib  --citeproc --csl bmc-bioinformatics.csl -f markdown out.md  --pdf-engine=xelatex -o paper.rtf
