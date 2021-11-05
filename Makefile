all: paper.pdf

refdb/ncbi_eukprot_met_arch_markers.fna:
	mkdir refdb
	cd refdb
	wget https://ndownloader.figshare.com/files/26173346
	tar -zxvf 26173346

refdbCrossValidation/nineTenth.fna.1.bt2: refdb/ncbi_eukprot_met_arch_markers.fna
	mkdir -pv refdbCrossValidation
	bash scripts/split_eukdb_reference.sh ./refdb/ ./refdbCrossValidation/oneTenth.fna ./refdbCrossValidation/nineTenth.fna
	bowtie2-build ./refdbCrossValidation/nineTenth.fna ./refdbCrossValidation/nineTenth.fna


tmp:
	mkdir -pv tmp

tmp/wgsimMutationRate.json: tmp refdb/ncbi_eukprot_met_arch_markers.fna
	python3 scripts/simulate_and_align.py --refdb-ncbi refdb/taxa.sqlite --refdb-marker-to-taxon-path refdb/busco_taxid_link.txt --reference refdb/ncbi_eukprot_met_arch_markers.fna --sim-source refdb/ncbi_eukprot_met_arch_markers.fna --dir tmp --verbose --read-lengths 100 --base-error-rates 0.0 --mutation-rates 0.0 0.001 0.01 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2 > tmp/wgsimMutationRate.json


tmpLeaveOneOut:
	mkdir -pv tmpLeaveOneOut

tmpLeaveOneOut/wgsimMutationRateLeaveOneOut.json: tmpLeaveOneOut refdbCrossValidation/nineTenth.fna.1.bt2
	python3 scripts/simulate_and_align.py --refdb-ncbi refdb/taxa.sqlite --refdb-marker-to-taxon-path refdb/busco_taxid_link.txt --reference refdbCrossValidation/nineTenth.fna --sim-source refdbCrossValidation/oneTenth.fna --dir tmpLeaveOneOut --verbose --read-lengths 100 --base-error-rates 0.0 --mutation-rates 0.0 0.01 0.025 0.05 > tmpLeaveOneOut/wgsimMutationRateLeaveOneOut.json

figures/wgsimMutationRate.png: tmp/wgsimMutationRate.json
	mkdir -pv figures
	python3 scripts/plot_recall_and_mapq_for_varying_mutation_rate.py --input-json tmp/wgsimMutationRate.json --output-png figures/wgsimMutationRate.png

figures/leaveOneOut.png: tmpLeaveOneOut/wgsimMutationRateLeaveOneOut.json
	mkdir -pv figures
	python3 scripts/plot_recall_and_mapq_for_varying_mutation_rate.py --input-json tmpLeaveOneOut/wgsimMutationRateLeaveOneOut.json --output-png figures/leaveOneOut.png

figures/bars.png: tmp/wgsimMutationRate.json
	mkdir -pv figures
	python3 scripts/plot_bars.py --input-json tmp/wgsimMutationRate.json --output-png figures/bars.png

figures/barsLeaveOneOut.png: tmpLeaveOneOut/wgsimMutationRateLeaveOneOut.json
	mkdir -pv figures
	python3 scripts/plot_bars.py --input-json tmpLeaveOneOut/wgsimMutationRateLeaveOneOut.json --output-png figures/barsLeaveOneOut.png 

supplement/wgsim.tsv: tmp/wgsimMutationRate.json
	mkdir -pv supplement
	python3 scripts/wgsim_to_tsv.py  --input-json tmp/wgsimMutationRate.json --output-tsv supplement/wgsim.tsv

supplement/wgsimLeaveOneOut.tsv: tmpLeaveOneOut/wgsimMutationRateLeaveOneOut.json
	mkdir -pv supplement
	python3 scripts/wgsim_to_tsv.py  --input-json tmpLeaveOneOut/wgsimMutationRateLeaveOneOut.json --output-tsv supplement/wgsimLeaveOneOut.tsv

paper.pdf: paper.md biblio.bib figures/wgsimMutationRate.png figures/leaveOneOut.png figures/bars.png figures/barsLeaveOneOut.png supplement/wgsim.tsv  supplement/wgsimLeaveOneOut.tsv
	pandoc -s --bibliography biblio.bib  --citeproc -f markdown paper.md -o paper.pdf
