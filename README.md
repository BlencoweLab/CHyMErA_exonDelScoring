# CHyMErA_exonDelScoring
Code to score CHyMErA exon deletions

Using read counts per guide pair as well as dropout scores from https://github.com/HenryWard/chymera-scoring, together with auxiliary information about exons such as expression and splicing status in the screened cell type, exons are scored for growth phenotypes based on the proportion of exon-deletion (intronic-intronic) guide pairs that drop out significantly.

**Input**:
- exonDeletionLibrary_normCounts_NovaSeq_18Sept18.txt.gz
