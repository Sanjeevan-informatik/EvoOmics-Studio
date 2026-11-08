# Roadmap

These items are planned, not delivered functionality.

1. **Genomics robustness:** htslib-backed BCF/VCF adapter; DP/GQ/FT thresholds; genomic site maps; indexed region access; haplotype-aware representation; sample metadata joins. Acceptance: fixture comparison against validated established readers and explicit allele/ploidy tests.
2. **Scale:** chunked Zarr/array storage, vectorized distances, memory/runtime benchmarks and resource estimates. Acceptance: representative barley datasets with measured peak memory and bounded chunk sizes.
3. **Phylogenetic validation:** neighbor joining, resampling support, established maximum-likelihood tool adapters and ascertainment-aware workflows. Acceptance: external reference results, fixed seeds and versioned tool logs. GTR+Gamma and likelihood scores will come from validated implementations, not placeholder formulas.
4. **Multi-omics:** separate count, methylation, abundance and sparse-matrix models, modality-specific QC and explicit sample linkage. Acceptance: known datasets and tests of normalization/missingness behavior.
5. **Workflow automation:** Snakemake/Nextflow templates and pinned environment specifications. Acceptance: a clean checkout reproduces documented outputs with no manual file edits.
6. **Studio UI:** project browser, QC, tree and matrix views backed by the existing service layer. Acceptance: usability checks plus protected job execution, cancellation, quotas and isolated data storage.

Before publication: choose a license, establish release tags/changelogs, add a software citation record with confirmed author details, benchmark realistic cohorts and conduct independent scientific review.
