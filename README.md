# Marker alignments paper 

This is a repository for a publication about CORRAL - [submitted to biorxiv](https://doi.org/10.1101/2022.03.09.483664). 

See also the repositories for:
- [`CORRAL`, the Nextflow pipeline](https://github.com/wbazant/CORRAL)
- [`marker_alignments`, a Python package](https://github.com/wbazant/marker_alignments)

It contains code needed to reproduce the analyses and simulations, as well as results.

## Main result files

### Simulated reads
Path: `supplement/simulatedReads.xlsx`

Contains a summary of our experiments with individual reads simulated, aligned to a reference, and tracked on a per-read level. There are two main variables: reads sampled from and aligned to original EukDetect reference or if read is sampled from hold-out set and aligned to a remaining set, and reads if are mutated before aligning. The results are summarized at dataset level for each run, and for the non-mutated reads we additionally provide a per-taxon breakdown.

### Simulated whole samples at low abundance
Path: `supplement/wgsimWholeSamplesLowAbundance.xlsx`

Contains a summary of results produced by our runs of CORRAL, EukDetect, and their variants, on a dataset of 335 simulated samples each with very low abundance.

### Simulated unknown species
Path: `supplement/wgsimWholeSamplesUnknownTaxaOneTenthCoverage.xlsx`

Contains a summary of results produced by our runs of CORRAL, EukDetect, and their variants, on a dataset of 338 simulated samples, at abundance corresponding to genome coverage of 0.1 but for species that are not in the reference.

### Simulated pairs of taxa
Path: `supplement/pairs.xlsx`

Contains a summary of results produced by our runs of CORRAL, EukDetect, and their variants, on simulated samples with pairs of species, in two subsets (a broad subset of 98 taxa and a focused subset of 50 taxa).

### DIABIMMUNE comparison
Path: `supplement/diabimmune.xlsx`

Contains our comparison of results on WGS taxa from the DIABIMMUNE study, between MicrobiomeDB and the original EukDetect publication

### MicrobiomeDB results
Path: `supplement/microbiomedb.xlsx`

Contains our summary of results from the release of MicrobiomeDB described in our publication.

## Code

This repository is also a Make pipeline. For any result here that interests you, you should be able to find in the Makefile the command for generating the file from scratch. You can also browse through the `scripts` directory.

### Subset selection
In a few places in our analysis we had to make choices about which taxa to analyse. You can find the procedures we used:

1. Splitting the reference into holdout and remaining sets is in `scripts/split_eukdb_reference.sh`
2. Selecting pairs of species for analysis is in `scripts/pick_pairs.py` 


## How to reproduce

Theoretically, you should be able to type one command:
```
make results
```
and re-create all our simulations and their analyses. This will take around 100 hours of compute on a standard laptop and about 100GB of temporary disk space. We did the analysis on a Ubuntu laptop, with some external dependencies including `bowtie2` and `wgsim`, Python, Perl and Bash.

If you attempt to do this and encounter any difficulties, please open an issue in this repository, we will be happy to help.

To recreate results for the DIABIMMUNE comparison or the MicrobiomeDB release, you need a distributed computing environment. To reproduce the results, run CORRAL with a list of samples as input.
