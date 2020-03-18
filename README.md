# CHyMErA_exonDelScoring
Code to score CHyMErA exon deletions

Using read counts per guide pair as well as dropout scores from https://github.com/HenryWard/chymera-scoring, together with auxiliary information about exons such as expression and splicing status in the screened cell type, exons are scored for growth phenotypes based on the proportion of exon-deletion (intronic-intronic) guide pairs that drop out significantly.

## Input
**Screen data**
- *exon_deletion_norm_counts.txt.gz*: Raw read counts and log2-fold changes per guide pair
- *intronic_intronic_guides.tsv.gz*: Mean log2-fold change and empirical FDR per intronic-intronic guide pair. See [here](https://github.com/HenryWard/chymera-scoring/blob/master/input/exon_deletion_norm_counts.txt) for derivation.
- *other_guides.tsv.gz*: Same as above, for other types of guide pairs

**Exon annotation**
- *Knockout.exons.human.csv.gz*: Exon annotation from screen design phase

**Additional information about gene expression and splicing in cell type**
This was obtained using [vast-tools](https://github.com/vastgroup/vast-tools) on RNA-Seq data from RPE1 cells, [GEO: GSE75189](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75189)
- *cRPKM_AND_COUNTS-Hsa1.tab.gz* Gene expression in RPE1
- *INCLUSION_LEVELS_FULL-Hsa1-hg19.tab.gz* Alternative splicing percent-spliced-in in RPE1

## Output
- *Fitness.RPE1.pdf*: Fitness effect (drop-out) of exonic Cas9 guides and annotation of essentiality cutoffs
- *Essentiality_RNAseq.pdf*: Mean fFitness effect of exonic Cas9 guides vs. expression of host genes and PSI of targeted alternative exons, showing that only expressed genes have fitness effects in the screen.
- *Knockout.exons_RPE1.pdf*: Numbers of exons by type according to an annotation modified by expression and splicing data from the same cell line.
- *ExonDeletionScores.RPE1.csv*: Table of signifcantly dropping out exon-targeting and control guide pairs, and hit annotation (Fig. 6a, left)
- *ScreenPerformance.RPE1.exonDelPairs.pdf*: Hit calling at various thresholds for proportion of guides that exhibit drop-out, for intronic-intronic guide pairs. (Fig. 6a, right)
- *ScreenPerformance.cutoffs.RPE1.singleCutPairs.pdf*: Hit calling at various thresholds for proportion of guides that exhibit drop-out, for intronic-intergenic control guide pairs.
- *HitCalling_0.18.pdf*: Rate of hits called at selected threshold of 18% of intronic-intronic guide pairs per exon (Fig. 6b)
- *Hits_RPE1_0.18.properties.pdf*: Fitness effect of exonic Cas9 guides for genes with significant exons (Fig. 6d) and exon legnth of significant and non-significant exons (Fig. S6a)
