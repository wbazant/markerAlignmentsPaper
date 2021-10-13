# Marker alignments paper
by Wojtek, Kathryn, Dan, possibly others

## Introduction

```
to do:
this is an in
put something here later
```


## Background
```
to do:
these are loose notes about of things I've heard about
research the topic
```

K-mer methods and de novo assembly are possible alternatives to read mapping.

Read mapping is the oldest, can be used with whole genomes.

Metaphlan is the canonical tool for running `bowtie2` on a reference of marker genes. Massive reference, [~1.1M unique clade-specific marker genes identified from ~100k reference genomes (~99,500 bacterial and archaeal and ~500 eukaryotic)](https://github.com/biobakery/biobakery/wiki/metaphlan3).

[Kaiju](https://kaiju.binf.ku.dk/) uses protein alignments, trying to get higher sensitivity.

EukDetect tries to be better at detecting just eukaryotes by lowering filtering criteria and carefully constructing a reference that skips the bacteria.


## Conclusion

```
to do:
this tries to explain why what I did is good and clever

use it to structure further work
```


Our contributions are:

- an explanation how off-target hits happen, which helps authors of tools that interpret read mapping as counts of taxa
- software that runs and interprets alignments to EukDetect's reference database, and possibly to other reference databases, which helps people who have metagenomic data and need to analyse in taxa present
- the analysis of MicrobiomeDB data, which helps people who want to know what eukaryotes are present in human microbiomes



We've not yet contributed, but could contribute:

- a reference for traces of animal DNA, which helps people who study samples where they might be present
- an analysis of host blood meals of a mosquito dataset, which helps epidemiologists and people who want to know what different mosquitoes feed on
- software for building references that will work well with our method, which helps bioinformaticians set up analyses like ours



We've had ideas about:

- using this work to detect anti-microbial resistance genes (is there a universal database of known AMR variants of genes?) 



Our work is new because:

- nobody else interprets low MAPQ reads as a good alignment in the wrong place
- nobody else interprets secondary alignments as evidence against taxon being present
- nobody else uses network methods


Our work is good because:

- the method works at very low abundances, it's even better at it than EukDetect
- the method reports mixtures of related species more sensitively than EukDetect
- the method does not skip, or bias counts against, species that only approximately match the reference

Our work would be even better if:
- the method reported estimates of the gap between reference and signal
- the method of reporting unknown taxa was better

## Methodology and results

### Definitions
```
to do:
this is a list of relevant words and how I use them
Remove the list, just introduce the words by correctly using them in sentences
```

**taxon** - an organism, in this context an organism that had its genome sequenced, appears in the reference database, and may or may not be present in the sequenced sample

**BUSCO** - a family of genes that are mostly present in each taxon and mostly single-copy

**marker** - a DNA sequence of a gene in a taxon that is assumed to be unique to the taxon. Markers in related taxa can be similar, for example if they belong to the same BUSCO

**reference database** - one of the inputs for an aligner, in this context it's a reference of markers that can be matched to

**alignment / hit** - a read in the sequenced sample that was found by an aligner to match a marker in the reference

Aligners like `bowtie2` can report multiple alignments per read. In that case, we can distinguish:

**primary alignment** - the best match for the read (based on alignment score, sequence identity, or other metric)

**secondary alignments** - alignments corresponding to matches for a read that are not the best one


We would also like to differentiate bewteen:

**target hit** - a read that comes from an organism $A$, and matches a marker $M'$ for a taxon $A'$, such that $A'$ is the closest taxon to $A$

**off-target hit** - a read that comes from an organism $A$, but matches a marker $M'$ for a taxon $A''$, and there is at least one taxon $A'$ closer to $A$ than $A''$ is



### Sources of off-target hits
```
to do:
this part is all anegdotal
it's about hits I saw in the data that I think are off-target and that I think I understand

add evidence to it or get rid of it
```

Even if the reference contains all taxa that might be in the sample, there are still some possible sources of off-target hits, like:

1. random bias: the sequencing process introduces errors and short sequences can coincide by chance

2. ubiquitous subsequences: in end-to-end alignments, these look like clipped alignments to where a reference finishes with a sequence corresponding to e.g. a binding site common in proteins

3. missing markers: reads for a marker that is missing in the reference instead align to similar markers

Our method addresses these kinds of off-target hits through a combination of filters or thresholds on read quality and length, alignment length, and numbers of markers required to detect a taxon.

An entirely different kind of off-target hits is due to a sequenced organism having no exact match in the reference. We will describe how this happens, and how it can be addressed by making use of secondary alignments.

#### A model for off-target hits due to inexact matches
```

to do:
this is my conjecture for how off-target hits can contain signal:
- an inexact match looks in alignments as a mix of related matches
- asking for best alignment sets up competitive mapping and that's not good
- secondary alignments are good for the next part


make it a cartoon or less long
merge with the next part (marker clusters)
show evidence that it happens, and that it happens in the way described
```

Suppose an organism $A$ has a version $b_A$ of a BUSCO $b$, and the reference contains markers $b_{A_1},  \dots  b_{A_n}$ for taxa $A_1,  \dots  A_n$. Let us say $A$ is most similar to $A1$ - perhaps it's another strain of the same species. Assume also a least common ancestor $A_0$ of $A$ and $A_1$, and $A_{00}$ of $A_0, A_2, \dots A_n$.
As mutations accumulate over time, we can predict $b_A$ will be most similar to $b_{A_1}$, but - in places where $A_1$ has diverged from $A_0$ - some segments of $b_A$ are most similar to other $b_{A_i}$. Some segments of $b_A$ could also be equally similar in all $b_{A_i}$, if there has been reason for that sequence not to change since the joint common ancestor $A_{00}$.

When reads from $b_A$ are aligned to each of the $b_{A_i}$, we expect match identity to form a distribution. The $A_i$ which differ more from $A$ should have to lower average match identity and fewer matches, but might still attract high identity matches.


So, competitive alignment of reads from $b_A$ between $b_{A_1} \dots b_{A_n}$ will not entirely favour the closest $A_1$.

If we ask an aligner to report a single best alignment for each read, we expect to see $h_1$ hits to $b_{A_1}$, and smaller amounts $h_i$ of hits for other $b_{A_i}$, such that $H = \sum_{i=1}^{n}h_i$ is proportional to the count of reads coming from sequencing b(A). We also expect ratios of $H$ to $h_i$ to be related to sequence similarity between $b_A$ and $b_{A_i}$ - the further $A$ is from $A_1$, the larger the number of off-target hits $H - h_1$.

The effect of $b_{A_1} \dots b_{A_n}$ 'competing' for the best alignment of each read is illuminated when an aligner is asked to report all reads. Some of $h_2  \dots  h_n$ are then accompanied by secondary alignments to $b_{A_1}$ - call them $s_1$, and some of $h_1$ will be accompanied by secondary alignments to $b_{A_2}  \dots  b_{A_n}$, $s_2 \dots s_n$. If $h_1$ is much larger than $h_i$, $s_1$ should be much smaller than $s_i$ and thus $\frac{h_1}{s_1}$ should be larger than $\frac{h_i}{s_i}$ and independent of $H$.

Thus secondary alignments help us differentiate the presence of an organism $A$ reported as many hits to $A_1$ and fewer hits to A2 from the presence of two unrelated organisms X and Y. For example, we can report a ratio of primary to secondary alignments for each taxon.

It is possible for $b_{A_1}$ to be very different from $b_A$, or missing from the reference entirely. Perhaps the genome of $A_1$ is incorrectly annotated, or $A_1$ has lost b when adapting to its niche.






### Building marker clusters
```
to do:
this introduces marker clusters

it tries to theorise what kinds of clusterings there can be
EukDetect uses a taxonomic clustering but we do an empirical one

remove most of it, it's incredibly wordy
add a visualisation of marker clusters in examples
```

The effect of including secondary alignments is to only add hits to sequences similar to ones already present - after all, they both match on a read. Because identifying off-target hits requires grouping similar markers, secondary alignments provide valuable context for what markers are generally similar to what is present in the sample.

In our method, we run an aligner with as many secondary alignments as we can computationally afford, and then build a similarity graph where all matched markers are nodes and counts of reads that align to both markers are weighted edges. We then pass the triples (marker1, marker2, weight) to a clustering program, MCL.

Using a machine learning program like MCL relieves us from more precise modelling of what it means for two markers to be similar, or relying on prior information on what should be matched together.

MCL produces clusters with several valuable properties:

- markers in a cluster generally come from a single BUSCO and closely related taxa, but not always 
- broadly similar markers are grouped in larger clusters
- unique hits correspond to clusters with single markers
- a marker sharing reads with multiple putative clusters either gets assigned to one of them, or results in two clusters being merged



### Reporting taxa
```
to do:
this is about how the software makes calls on the taxon level
currently we sequentially apply a few filters and transforms - this doesn't even describe them all
the unknown taxon transform doesn't work as well as it should yet

get some clarity on the method and try to improve it?
```

#### Filter built from marker clusters
Clusterings of markers which attract similar matches allows for classifying taxa by having each marker cluster "vote" for its taxon, based on how that taxon's markers do in that cluster. This can be done in a number of ways as long as the marker clusters are required to only be approximately correct.

We choose to measure how good an alignment is through match identity - a number of bases that agree between query and reference divided by the alignment length. Then for each marker cluster, we compute an average for all matches in the cluster, as well as an average for matches in each marker. Then we discard taxa where less than half of each taxon's markers that are at least average in their cluster.

We can expect that this filter will work well enough to discard the additional taxa introduced to the result when setting the aligner to report secondary alignments, since they are on average inferior. Similarly, a version of a BUSCO that is overall inferior but has a locally better subsequence can be expected to accrue bad matches, get paired up with overall better versions of the BUSCO in the marker clustering process, and then help get its taxon rejected.

The filter based on marker clusters leaves us with taxa whose markers are not mostly off-targets of other markers in the reference. There can still be a gap between reference and signal we might get a mixture of matches at a marker level.

#### Taxon transform + identity filter
We build another graph where taxa are nodes and (directed) edges is a fraction of shared reads. Then we run MCL again to get taxon clusters.

This groups taxa which are similar and, hopefully, each corresponds to at most one taxon. 

#### Minimal evidence filter

After that there's still a need for general filtering like EukDetect and Metaphlan do.

We do it at the end and have separate thresholds for called taxa and vaguely identified blobs, as the burden of evidence is higher for taxa with lower identity.

This part is a combination of read counts, marker counts.

### Implementation
```
to do:
This section is about the software
Figure out what should actually be here
```

We implement our way of interpreting results of alignments to markers as a Python package, `marker_alignments`. It uses a module `pysam` to read alignments into a following SQLite table:

| column | type | description |
| -- | -- | -- |
| query | text | read identifier |
| taxon | text | name of matched taxon |
| marker | text | name of matched taxon |
| coverage | number | fraction of marker covered by match  | 
| identity | number | fraction of bases agreeing between query and reference  | 

Counts of entries and the coverage field are used for quantification, and identity is mainly used for clustering and filtering.

We implement filtering and reporting using SQL queries and Python code. Additionally, clustering uses additional tables and a `markov_clustering` package, a Python implementation of the MCL algorithm.

In addition we also provide a Nextflow workflow, `marker-alignments-nextflow`.

### Shortcomings and potential future work
```
to do:
This section is a catch-all of stuff that is bad
replace it with an entirely different section on the same topic
```

We add an up-front threshold on alignment length of 60 bases. This is a filter borrowed from EukDetect and we do not have a good rationale for keeping the filter, except without the filter, some samples (our analysis of NICU NEC dataset) contain a lot of matches to taxa that are not plausibly present given the sequenced samples. We suspect they might be incorrect annotations, or annotations that end with common elements (like binding sites) but we have not investigated further.

We also see an occasional inclusion of matches to model organisms. We suspect these are BUSCOs that are commonly present, but not commonly annotated. They happen rarely enough to occur twice for any dataset, so an ad hoc criterion of (TODO implement this!) twenty reads for taxa reported only in one sample per dataset removes these matches and keeps matches where we have much more confidence despite them being present only in one dataset.

## Examples
```
to do:
this is a catch-all with data examples
sprinkle these throughout, remove the rest
```

### Example 1 - off-target hits in more than just a genus
SRR6262267 is a run from a sample dominated by *Trichosporon asahii* - according to SRA Traces, 23.68% of the reads in the sample can be attributed to this organism ([source](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR6262267)).

```
TODO an example - perhaps a heatmap of counts for BUSCOs in each Trichosporon?
Actually, EukDetect reports only T. asahii, because the other off-target hits are for the same genus - find something :).
TODO this will require some visualisation tools.
```

### Example 2 - sticking to the reference too closely brings up nothing
Mucor example, demonstrate EukDetect returns nothing which it really should.

### Example 3 - BUSCOs vs clusters
Maybe a drawing or a visualisation: a graph with BUSCOs corresponding to a shading of each node, added edges, and a shading or a line around clusters that end up together?

## Applications

### MicrobiomeDB studies
```
to do:
This is about MicrobiomeDB results.
Decide whether it should even be here.
```

We analyzed all data on MicrobiomeDB.

| study | num samples | num samples with reported taxa | num reported taxa |
| -- | -- |-- | -- |
| HMP | xxx | yyy | zzz |

The data can be viewed on MicrobiomeDB :).

#### DIABIMMUNE
```
to do:
This is an expanded part of the above that does a comparison for DIABIMMUNE.
Write a program to report the numbers, for transparency and because you will need to redo it.
```

Here's a not yet formal comparison to DIABIMMUNE.

Just the paper DIABIMMUNE reported 8 samples. I think it's mostly the difference between references - new one has a bit more.
We're more courageous in reporting species, and report 119 extra results. It shows up for the most common S. cerevisiae - we call it 42 more times - and M. restricta - 9 more times.

#### 3100266
```
> EukDetect reports P. nordicum, but we realise it's a cluster:
id      taxon   num_markers     num_reads       avg_identity
1       1429867|Penicillium_camemberti_FM_013   34      73      0.960183806201549
1       1439350|Penicillium_fuscoglaucum_FM041  33      74      0.962397902985074
1       1931374|Penicillium_sp_BW_MB    101     243     0.958158950450451
1       1931375|Penicillium_sp_BW_12    10      31      0.959146928571428
1       229535|Penicillium_nordicum     45      102     0.954239219101123
1       2488753|Penicillium_sp_SPG-F1   10      29      0.964647588235294
1       2593313|Penicillium     124     313     0.967739331578951
1       48697|Penicillium_freii 117     288     0.959325344168261
1       5073|Penicillium        306     359     0.962057904903424
1       60169|Penicillium_polonicum     41      78      0.959593049295774
1       60171|Penicillium_verrucosum    40      70      0.955100714285714
1       60172|Penicillium_solitum       39      81      0.959004702702702
2       5480|Candida_parapsilosis       8       10      0.996333894736842
```
and because they're all above 95%, we return them all. What I would prefer is to return the whole cluster, but how should I know?

#### 3104340
```
G80329
diabimmune-paper.tsv.csv:G80329,Candida parapsilosis,5480,27,78,15.70%,8.04%,99.64%
diabimmune-paper.tsv.csv:G80329,Rhodotorula,5533,3,5,3.53%,3.76%,100.00%
diabimmune-paper.tsv.csv:G80329,Chaetomium globosum CBS 148.51,306901,3,5,1.69%,4.14%,99.79%
diabimmune.tsv:3104340  G80329  237561  Candida_albicans_SC5314 8
diabimmune.tsv:3104340  G80329  306901  Chaetomium_globosum_CBS_14851   3
diabimmune.tsv:3104340  G80329  5480    Candida_parapsilosis    19
diabimmune.tsv:3104340  G80329  5533    Rhodotorula     3
```
we are able to show that the presence of additional C. albicans is very convincing - the hits are in entirely different marker clusters. (edited) 


### Potential application 1 - host blood meal
Identify host blood meal in metagenomic studies of mosquitos.


### Potential application 2 - AMR
Detection/quantification of antimicrobial resistance (AMR) genes in metagenomic samples

\bibliography{paper}

