all: paper.pdf

refdb/ncbi_eukprot_met_arch_markers.fna:
	mkdir refdb
	cd refdb
	wget https://ndownloader.figshare.com/files/26173346
	tar -zxvf 26173346

tmp:
	mkdir -pv tmp

tmp/wgsimMutationRate.json: tmp refdb/ncbi_eukprot_met_arch_markers.fna
	python3 scripts/simulate_and_align.py --reference refdb/ncbi_eukprot_met_arch_markers.fna --sim-source refdb/ncbi_eukprot_met_arch_markers.fna --dir tmp --verbose --read-lengths 100 --base-error-rates 0.0 --mutation-rates 0.0 0.001 0.01 0.05 0.075 0.1 0.125 0.15 0.2 > tmp/wgsimMutationRate.json

figures:
	mkdir -pv figures

figures/wgsimMutationRate.png: figures tmp/wgsimMutationRate.json
	python3 scripts/plot_recall_and_mapq_for_varying_mutation_rate.py --input-json tmp/wgsimMutationRate.json --output-png figures/wgsimMutationRate.png

paper.pdf: paper.md biblio.bib figures/wgsimMutationRate.png
	pandoc -s --bibliography biblio.bib  --citeproc -f markdown paper.md -o paper.pdf
