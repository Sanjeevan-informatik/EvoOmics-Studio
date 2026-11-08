# Contributing

Install with `python -m pip install -e .` and run `python -m unittest discover -s tests -v`. Keep parsing, scientific algorithms and orchestration separate. Add analytical or externally validated expected results for algorithm changes. Report missingness, units, coordinate conventions and assumptions explicitly. Do not commit private datasets, credentials, generated build bundles or run outputs. Synthetic fixtures should stay small and explain the scientific edge case they cover.

Changes to output structure require a manifest schema review. New features must update the capability table and distinguish working behavior from planned integration. Discuss licensing with the repository owner before importing third-party code.
