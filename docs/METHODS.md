# Scientific methods

For each sample pair, use only positions where both characters belong to A/C/G/T. Let n be that overlap, d the mismatch count, t the transition count and v the transversion count.

- p-distance = d / n.
- JC69 = -3/4 × ln(1 - 4p/3), only when p < 3/4.
- K2P = -1/2 × ln(1 - 2P - Q) - 1/4 × ln(1 - 2Q), with P=t/n and Q=v/n. Both logarithm arguments must be positive.

Transitions are A↔G and C↔T. K2P distinguishes transitions and transversions but assumes equal equilibrium base frequencies. JC69 assumes equal base frequencies and equal substitution rates. Neither model corrects SNP ascertainment in this implementation.

UPGMA repeatedly merges the closest clusters. The new cluster's distance to another cluster is the average weighted by the number of original samples in each merged cluster. Node height is half the merge distance; a child branch length is parent height minus child height. Zero-length branches are valid. Deterministic tie handling depends on input order and is not biological evidence.

Example: distances AB=2, AC=4, BC=4 merge A/B at height 1 and join C at height 2. This is tested exactly. JC69 at p=0.25 and K2P at P=Q=0.1 are checked against independent analytical expressions in the tests.

A molecular clock is assumed by UPGMA, not established by using it. No support values are estimated. Differences in missing-data patterns, ascertainment, linkage, recombination and evolutionary rates may make the resulting tree misleading. Use the exported QC and distances to judge whether an analysis is appropriate, and validate research claims with established phylogenetic tools and appropriate site ascertainment treatment.
