# Data contracts — version 1.0

## Implemented DNA alignment

At least two unique sample IDs, restricted to ASCII letters, numbers, dot, underscore and hyphen. FASTA uses the first whitespace-delimited header token as the ID. Sequences have equal positive length; uppercase DNA and IUPAC ambiguity symbols, N, ? and - are accepted (FASTA is normalized to uppercase). Only A/C/G/T are callable for distances.

VCF input uses the standard nine fixed columns followed by at least two samples. For retained PASS biallelic SNPs, GT must be haploid or diploid with allele indices 0/1 or missing dots; FORMAT may place GT at any index. A missing GT subfield is masked. Duplicate retained chromosome/position pairs fail. VCF loci are represented in input order; `retained_columns.json` maps back to the retained-SNP import order, not all original VCF lines. Recover genomic coordinates from the hashed original VCF using the documented import rules. No site sorting, annotation or genome-build conversion is performed.

The VCF reader intentionally implements a constrained subset, not full VCF specification validation. Upstream validation and caller-specific QC remain necessary.

## Run manifest

`schema_version`, software/version, Python version, UTC creation time, input filename/hash/format, all analysis parameters, artifact filenames/hashes and warnings. Relative artifact paths refer to the run directory. The manifest cannot hash itself. Timestamps vary across runs; scientific result files are deterministic for identical input order and parameters. Provenance includes neither a full environment lock nor a source commit hash; use a released source revision alongside the manifest for archival studies.

## Planned modality contracts (not supported yet)

| Modality | Native values | Required future metadata/QC |
|---|---|---|
| RNA-seq | Integer counts or explicitly labeled TPM | Gene ID namespace, library size, normalization, batch |
| Methylation | Beta in [0,1] or declared M-values | Probe/locus identifiers, detection QC, platform, genome build |
| Proteomics | Instrument intensity or declared normalized abundance | Protein IDs, missingness, batch, transformation |
| Single-cell | Sparse counts | Cell/sample mapping, feature IDs, doublet and library QC |
| Annotation | Genomic features | Assembly, coordinate convention, annotation release |

Integration will require explicit sample joins, unmatched-ID reports and appropriate normalization. No automatic equivalence between measurement units is assumed.
