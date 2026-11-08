# Architecture

## Boundaries

`io` parses external representations into the validated `Alignment` object. `qc` applies explicit site-callability rules. `phylogeny` performs distance calculation and tree inference without filesystem operations. `pipeline.run` coordinates these components and writes auditable outputs. `cli` handles arguments and friendly errors. A future UI or workflow runner should call the same service, not reimplement scientific logic.

| Location | Responsibility |
|---|---|
| `src/evoomics/models.py` | Domain types and validation errors |
| `src/evoomics/io/` | FASTA/VCF input adapters |
| `src/evoomics/qc.py` | Site filtering and sample QC |
| `src/evoomics/phylogeny/` | Pure distance and UPGMA algorithms |
| `src/evoomics/pipeline.py` | Run lifecycle, export and provenance |
| `src/evoomics/cli.py` | User command interface |
| `examples/` | Small synthetic fixtures, no personal or clinical data |
| `tests/` | Scientific unit tests and end-to-end regressions |
| `.github/workflows/` | Automated package and regression checks |

## Run lifecycle

1. Validate requested parameters and reject existing output directories.
2. Hash input, parse records, validate alignment and apply site QC.
3. Compute distances and comparable-site counts; stop on insufficient overlap or saturation.
4. Build UPGMA and check that the input did not change while being read.
5. Write results in a temporary staging directory and hash every result.
6. Exclusively create the requested directory and move files into it, cleaning up on ordinary errors.

Output publication is not a database transaction: an abrupt process kill during the final move can leave an incomplete run. Consumers must require `provenance.json` and verify every listed hash. Files are local; there is no upload endpoint, authentication system, network job queue or remote command execution.

## Extension design

Future adapters must preserve sample identifiers and keep modality-specific units and QC. Do not coerce count matrices, methylation beta values or protein intensities into nucleotide alignments. Introduce separate domain types when a supported use case is implemented. Planned analysis modules need an explicit input contract, versioned parameters, scientific validation, failure semantics and resource limits before being marked available.

For scale, retain the CLI/service boundary while replacing the in-memory alignment with chunked arrays and adding background jobs. For external tools, use argument lists, pinned tool versions, bounded resource usage and captured logs; never concatenate user input into shell commands. A public API requires authentication, quotas, isolated storage and workers before exposing compute jobs.
