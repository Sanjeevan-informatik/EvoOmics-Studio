# EvoOmics Studio

A reproducible foundation for evolutionary genomics, with a modular path toward multi-omics research.

**Version 0.1.0: working Python library and CLI.** Start with aligned DNA FASTA or a filtered SNP VCF and obtain a distance matrix, exploratory UPGMA tree, QC summaries and a provenance manifest. No web interface or multi-omics inference is implemented yet. Research software; not validated for clinical decisions.

## Quick start

Python 3.11 or newer. The analysis package has **no third-party runtime dependencies**.

```bash
git clone https://github.com/Sanjeevan-informatik/EvoOmics-Studio.git
cd EvoOmics-Studio
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
evoomics analyze --input examples/aligned.fasta --format fasta --model jc69 --output results/fasta-demo
evoomics analyze --input examples/tiny.vcf --format vcf --output results/vcf-demo
```

Windows PowerShell: replace the activation command with `.\.venv\Scripts\Activate.ps1`.

Offline, without installing the package (Bash):

```bash
PYTHONPATH=src python -m evoomics analyze --input examples/aligned.fasta --format fasta --output results/offline-demo
```

Use a **new output directory for each run**. Existing directories are never overwritten.

## What works today

| Capability | Behavior |
|---|---|
| Aligned DNA FASTA / FASTA.gz | Validates unique IDs, alphabet and equal nonzero lengths |
| VCF / VCF.gz | PASS biallelic SNPs; haploid and diploid GT; phased or unphased input |
| Heterozygotes | Mask to N by default; optional IUPAC export; both excluded from distance calculations |
| Missingness QC | Per-sample callability, site filtering and comparable-site counts |
| Distances | p-distance, JC69 and K2P with saturation checks |
| Trees | Deterministic, cluster-size weighted UPGMA; Newick output |
| Reproducibility | Input/output SHA-256, parameters, software/Python versions, UTC timestamp |
| Engineering | Separated adapters, models, QC, scientific algorithms, pipeline and CLI; tests and CI |

Example with stricter site/overlap requirements (choose thresholds for your dataset):

```bash
evoomics analyze --input examples/aligned.fasta --format fasta --output results/strict --max-site-missing 0 --min-comparable-sites 10
```

## Outputs

| File | Meaning |
|---|---|
| `alignment.fasta` | Imported DNA representation after site filtering |
| `distances.tsv` | Symmetric pairwise distances; p is a proportion, JC69/K2P estimate substitutions/site under their assumptions |
| `comparable_sites.tsv` | Callable A/C/G/T overlap for every sample pair |
| `tree.nwk` | UPGMA dendrogram, with branch lengths in distance units |
| `qc.json` | Input filtering counts and post-filter sample callability |
| `retained_columns.json` | Zero-based column indices in the imported alignment; not genomic coordinates |
| `provenance.json` | Versioned manifest, hashes, parameters and warnings |
| `report.md` | Short run summary and scientific limitations |

## Scientific boundaries

- Input FASTA **must already be aligned**. Equal length is necessary but does not establish homology. RNA `U` and protein sequences are rejected.
- Missing, gapped and ambiguous bases are excluded using **pairwise deletion**. Differing denominators can distort comparisons; inspect overlap counts and consider complete-site filtering.
- A VCF consensus is not a reconstructed haplotype. Phasing is not used; heterozygosity is not converted to a random allele. Multiallelic records, indels and non-PASS records (including `.`) are skipped and counted. Invalid retained SNP genotypes fail explicitly.
- DP/GQ and per-sample FT thresholds are **not** applied. Perform appropriate upstream variant/genotype QC before import. BCF, gVCF blocks, polyploid GT and reference reconstruction are unsupported.
- SNP-only data omit invariant sites and cause **ascertainment bias**. Corrected distances on VCF are blocked unless `--allow-snp-correction` is explicitly supplied. That flag acknowledges bias; it does not fix it. SNP-only FASTA has the same limitation, which cannot be inferred from its file format.
- UPGMA assumes equal evolutionary rates (a molecular clock). This is an exploratory tree, without bootstrap support, likelihood scores, model selection or a validated evolutionary root. Ties are resolved deterministically by input order.
- Current limits: **200 samples and 100,000 imported sites**. These are guardrails, not performance guarantees. Data are held in memory; distance computation is O(samples² × sites), and the simple UPGMA implementation is O(samples³). Large barley cohorts need chunked storage and optimized algorithms.

See [scientific methods](docs/METHODS.md), [architecture](docs/ARCHITECTURE.md), [data contracts](docs/DATA_CONTRACTS.md) and [roadmap](docs/ROADMAP.md).

## Development

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
```

Tests cover analytical distances, saturation, missingness, genotype interpretation, weighted UPGMA, reproducibility, manifest checksums, compressed input and CLI error handling. GitHub Actions runs the suite on Python 3.11–3.13 and verifies a wheel installation in a clean virtual environment.

Optional container:

```bash
docker build -t evoomics-studio .
docker run --rm -v "$PWD/examples:/data:ro" -v "$PWD/results:/results" evoomics-studio analyze --input /data/aligned.fasta --format fasta --output /results/container-demo
```

The container runs as UID 10001. Prepare a writable output mount for that user. The container has not been benchmarked or hardened as a multi-user service.

## Project history

The supplied VCFToPhylogenetic archive informed the scope: genotype-based tree inference and reproducible phylogenetic analysis. Its README acknowledges IPK DivBrowse and supervision by Patrick Koenig. This foundation is newly written; legacy code, compiled bytecode and generated frontend bundles were not imported. No license was present in the supplied archive. A license for this repository should be selected by the owner before third-party redistribution.
