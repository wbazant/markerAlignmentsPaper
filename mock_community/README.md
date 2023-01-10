# How to reproduce
Most of the work for reproducing the publication can be done on a laptop. This part involved a server, but still no distributed computing environment is needed.

### Download .fastqs
We downloaded fastqs from PRJEB38036 that correspond to mock community samples using `wget` from ENA.

The compressed .fastqs took 344G of disk space, and took about a week, at speed just under 0.7MBps.

### Download EukDetect's database

### Run `bowtie2` twice
We run it once in best-alignment mode, and once to report multiple alignments for CORRAL (specified as -k10 or top ten).

### Fetch .sams from the server
Zipped and fetched with rsync

### Run the `marker_alignments` wrapper
We wrote and then used the `do_all_marker_alignments.sh` wrapper available in this folder.

### Patch EukDetect to allow alignment files as input
We replaced a part of the pipeline doing the alignments, with a `cat {input.sam}` using the reaady alignments. See `eukdetect-bowtie2-input.rules` and `conf.yaml`.

### Run patched EukDetect
This was the command:
```
snakemake --snakefile ./eukdetect-bowtie2-input.rules --configfile ./conf.yaml  --cores 1 runall
```

### Make auxiliary files
The results reference the original publication of the data. Go to [Yang et al. (2020)](https://doi.org/10.1093/gigascience/giaa071) , download the `Supplementary Tables S1-9.xlsx`, and unpack here. The other auxiliary files are provided for convenience: `errs.tsv`, `alignment_counts_all.tsv`, `alignment_counts_best.tsv`, `read_counts.tsv`.

### Create the spreadsheet
Command:
```
python3 make_spreadsheet.py
```


