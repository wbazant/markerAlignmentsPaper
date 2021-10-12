\title{Marker alignments paper}
\author{You}

## Introduction

## Prior work 
EukDetect \cite{lind2021accurate}

## About the method

### Definitions
**taxon** - an organism, in this context an organism that had its genome sequenced, appears in the reference database, and may or may not be present in the sequenced sample

**BUSCO** - a family of genes that are mostly present in each taxon and mostly single-copy

**marker** - a DNA sequence of a gene in a taxon that is assumed to be unique to the taxon. Markers in related taxa can be similar, for example if they belong to the same BUSCO

**reference database** - one of the inputs for an aligner, in this context it's a reference of markers that can be matched to

**alignment / hit** - a read in the sequenced sample that was found by an aligner to match a marker in the reference

Aligners like `bowtie2` can report multiple alignments per read. In that case, we can distinguish:

**primary alignment** - the best match for the read (based on alignment score, sequence identity, or other metric)

**secondary alignments** - alignments corresponding to matches for a read that are not the best one


We would also like to differentiate bewteen:

**target hit** - a read that comes from an orgamism $A$, and matches a marker $M'$ for a taxon $A'$, such that $A'$ is the closest taxon to $A$

**off-target hit** - a read that comes from an organism $A$, but matches a marker $M'$ for a taxon $A''$, and there is at least one taxon $A'$ closer to $A$ than $A''$ is



### Sources of off-target hits
Off-target hits could happen due to random bias - the sequencing process introduces errors and short sequences can coincide by chance.

#### A model for off-target hits
Suppose an organism $A$ has a version $b_A$ of a BUSCO $b$, and the reference contains markers $b_{A_1},  \dots  b_{A_n}$ for taxa $A_1,  \dots  A_n$. Let us say $A$ is most similar to $A1$ - perhaps it's another strain of the same species. Assume also a least common ancestor $A_0$ of $A$ and $A_1$, and $A_{00}$ of $A_0, A_2, \dots A_n$.
As mutations accumulate over time, we can predict $b_A$ will be most similar to $b_{A_1}$, but - in places where $A_1$ has diverged from $A_0$ - some segments of $b_A$ are most similar to other $b_{A_i}$. Some segments of $b_A$ could also be equally similar in all $b_{A_i}$, if there has been reason for that sequence not to change since the joint common ancestor $A_{00}$.

This generates a prediction of what to expect when reads from $b_A$ are aligned to each of the $b_{A_i}$ - match identity should form a distribution, with increased evolutionary distance between $A$ and $A_i$ leading to lower average identity and fewer matches. Additionally, competitive alignment of reads from $b_A$ between $b_{A_1} \dots b_{A_n}$ will not entirely favour the closest $A_1$.

If we ask an aligner to report a single best alignment for each read, we expect to see $h_1$ hits to $b_{A_1}$, and smaller amounts $h_i$ of hits for other $b_{A_i}$, such that $H = \sum_{i=1}^{n}h_i$ is proportional to the count of reads coming from sequencing b(A). We also expect ratios of $H$ to $h_i$ to be related to sequence similarity between $b_A$ and $b_{A_i}$ - the further $A$ is from $A_1$, the larger the number of off-target hits $H - h_1$.

The effect of $b_{A_1} \dots b_{A_n}$ 'competing' for the best alignment of each read is illuminated when an aligner is asked to report all reads. Some of $h_2  \dots  h_n$ are then accompanied by secondary alignments to $b_{A_1}$ - call them $s_1$, and some of $h_1$ will be accompanied by secondary alignments to $b_{A_2}  \dots  b_{A_n}$, $s_2 \dots s_n$. If $h_1$ is much larger than $h_i$, $s_1$ should be much smaller than $s_i$ and thus $\frac{h_1}{s_1}$ should be larger than $\frac{h_i}{s_i}$ and independent of $H$.

Thus secondary alignments help us differentiate the presence of an organism $A$ reported as many hits to $A_1$ and fewer hits to A2 from the presence of two unrelated organisms X and Y. For example, we can report a ratio of primary to secondary alignments for each taxon.

It is possible for $b_{A_1}$ to be very different from $b_A$, or missing from the reference entirely. Perhaps the genome of $A_1$ is incorrectly annotated, or $A_1$ has lost b when adapting to its niche.


#### Example 1 - off-target hits in more than just a genus
SRR6262267 is a run from a sample dominated by *Trichosporon asahii* - according to SRA Traces, 23.68% of the reads in the sample can be attributed to this organism ([source](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR6262267)).

TODO an example - perhaps a heatmap of counts for BUSCOs in each Trichosporon?
Actually, EukDetect reports only T. asahii, because the other off-target hits are for the same genus - find something :).
TODO this will require some visualisation tools.

#### Example 2 - sticking to the reference too closely brings up nothing
Mucor example, demonstrate EukDetect returns nothing which it really should.


### A case for using all alignments
In the model described above, off-target hits happen when an organism in the sample has a version of a BUSCO which is unexactly similar to multiple markers, and this appears in the alignments as a mix of hits. Removing off-target hits requires a grouping of markers. EukDetect uses filtering on alignment length and the MAPQ field (which skews the results further away from off-target hits) followed by a taxonomic grouping of taxa by genus and a comparison on sequence identity.

This approach has a potential cost of sacrificing signal in clipped alignments at the ends of marker sequences (where a read might contain sequence past the end of the marker), and in hits in regions where markers share a sequence (The MAPQ field is defined by SAM manual as an expression of probability that the mapping position is wrong). Additionally, rejecting off-target hits has to be reconciled with sensitivity to presence of multiple related species in the sample.

A method of summarising alignments that can make use of alignments with worse match identity can extract more information from the sequencing data. If we can treat a match with lower identity is evidence for presence of something other than the matched taxon, the more alignments the better - and we only would only need to balance this benefit against its dimishing returns and rising computational cost to obtaining more alignments.

### Marker clusters
The effect of including secondary alignments is to only add hits to sequences similar to ones already present - after all, they both match on a read. Because identifying off-target hits requires grouping similar markers, secondary alignments provide valuable context for what markers are generally similar to what is present in the sample.

In our method, we run an aligner with as many secondary alignments as we can computationally afford, and then build a similarity graph where all matched markers are nodes and counts of reads that align to both markers are weighted edges. We then pass the triples (marker1, marker2, weight) to a clustering program, MCL.

Using a machine learning program like MCL relieves us from more precise modelling of what it means for two markers to be similar, or relying on prior information on what should be matched together.

MCL produces clusters with several valuable properties:
- markers in a cluster generally come from a single BUSCO and closely related taxa, but not always 
- broadly similar markers are grouped in larger clusters
- unique hits correspond to clusters with single markers
- a marker sharing reads with multiple putative clusters either gets assigned to one of them, or results in two clusters being merged


#### Example 3
TODO maybe a drawing or a visualisation: a graph with BUSCOs corresponding to a shading of each node, added edges, and a shading or a line around clusters that end up together?

### Using marker clusters to identify off-target hits
Clusterings of markers which attract similar matches allows for classifying taxa by having each marker cluster "vote" for its taxon, based on how that taxon's markers do in that cluster. This can be done in a number of ways as long as the marker clusters are required to only be approximately correct.

We choose to measure how good an alignment is through match identity - a number of bases that agree between query and reference divided by the alignment length. Then for each marker cluster, we compute an average for all matches in the cluster, as well as an average for matches in each marker. Then we discard taxa where less than half of each taxon's markers that are at least average in their cluster.

We can expect that this filter will work well enough to discard the additional taxa introduced to the result when setting the aligner to report secondary alignments, since they are on average inferior. Similarly, a version of a BUSCO that is overall inferior but has a locally better subsequence can be expected to accrue bad matches, get paired up with overall better versions of the BUSCO in the marker clustering process, and then help get its taxon rejected.

TODO a diagram or example

### Implementation

Here we pro

### Comparison with EukDetect
We apply our

### TODO summary of data loaded? Or maybe one good dataset
- how many samples, whi

### Shortcomings and potential future work
We add an up-front threshold on alignment length of 60 bases. This is a filter borrowed from EukDetect and we do not have a good rationale for keeping the filter, except without the filter, some samples (our analysis of NICU NEC dataset) contain a lot of matches to taxa that are not plausibly present given the sequenced samples. We suspect they might be incorrect annotations, or annotations that end with common elements (like binding sites) but we have not investigated further.

We also see an occasional inclusion of matches to model organisms. We suspect these are BUSCOs that are commonly present, but not commonly annotated. They happen rarely enough to occur twice for any dataset, so an ad hoc criterion of (TODO implement this!) twenty reads for taxa reported only in one sample per dataset removes these matches and keeps matches where we have much more confidence despite them being present only in one dataset.



\bibliography{paper}

