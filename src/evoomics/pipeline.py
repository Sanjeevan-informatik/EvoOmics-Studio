"""Reusable analysis service, shared by the CLI and future workflow adapters."""
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import shutil
import tempfile

from evoomics import __version__
from evoomics.io.sequences import read_fasta, read_vcf
from evoomics.models import ValidationError
from evoomics.phylogeny.distance import MODELS, distance_matrix
from evoomics.phylogeny.tree import upgma
from evoomics.qc import filter_alignment


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def write_matrix(path, samples, matrix):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id", *samples])
        for name, values in zip(samples, matrix):
            writer.writerow([name, *values])


def run(input_path: Path, output_path: Path, *, input_format="fasta", model="p",
        heterozygotes="mask", max_site_missing=1.0, min_comparable_sites=1,
        allow_snp_correction=False):
    input_path, output_path = Path(input_path), Path(output_path)
    if input_format not in {"fasta", "vcf"} or model not in MODELS:
        raise ValidationError("Unknown input format or distance model")
    if min_comparable_sites < 1:
        raise ValidationError("min_comparable_sites must be positive")
    if output_path.exists():
        raise ValidationError("Output already exists; choose a new run directory")
    if input_format == "vcf" and model != "p" and not allow_snp_correction:
        raise ValidationError("SNP-only VCF has ascertainment bias; use p or explicitly --allow-snp-correction")
    before = sha256(input_path)
    alignment, ingestion = (read_fasta(input_path) if input_format == "fasta"
                            else read_vcf(input_path, heterozygotes))
    filtered, qc, keep = filter_alignment(alignment, max_site_missing)
    matrix, counts = distance_matrix(filtered, model, min_comparable_sites)
    newick = upgma(filtered.samples, matrix)
    if sha256(input_path) != before:
        raise ValidationError("Input changed during analysis; retry with an immutable input")
    warnings = [
        "Exploratory UPGMA assumes a molecular clock; no bootstrap support or likelihood inference.",
        "Pairwise deletion may use different sites per pair; inspect comparable_sites.tsv.",
        "FASTA must be aligned; equal sequence lengths do not prove biological homology.",
    ]
    if input_format == "vcf":
        warnings += ["VCF output contains retained variant records only, not a whole-genome alignment.",
                     "Distances on SNP-ascertained data are not genome-wide substitution estimates.",
                     "Only PASS records are accepted; DP/GQ and sample FT filters are not applied."]
    else:
        warnings += ["SNP-only FASTA has the same ascertainment bias as SNP-only VCF."]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".evoomics-", dir=output_path.parent))
    try:
        (stage / "alignment.fasta").write_text("".join(
            f">{name}\n{seq}\n" for name, seq in zip(filtered.samples, filtered.sequences)), encoding="utf-8")
        write_matrix(stage / "distances.tsv", filtered.samples, matrix)
        write_matrix(stage / "comparable_sites.tsv", filtered.samples, counts)
        (stage / "tree.nwk").write_text(newick + "\n", encoding="utf-8")
        write_json(stage / "qc.json", {"ingestion": ingestion, "alignment": qc})
        # Positions refer to the imported alignment, not genomic VCF coordinates.
        write_json(stage / "retained_columns.json", {"index_basis": "zero-based imported alignment", "columns": keep})
        (stage / "report.md").write_text(
            f"# EvoOmics Studio analysis\n\nSamples: {len(filtered.samples)}\n\n"
            f"Retained sites: {filtered.sites}\n\nDistance: {model}; tree: UPGMA\n\n"
            + "\n".join(f"- {w}" for w in warnings) + "\n", encoding="utf-8")
        artifacts = {p.name: sha256(p) for p in sorted(stage.iterdir())}
        write_json(stage / "provenance.json", {
            "schema_version": "1.0", "software": "evoomics-studio", "version": __version__,
            "python": platform.python_version(), "created_utc": datetime.now(timezone.utc).isoformat(),
            "input": {"filename": input_path.name, "sha256": before, "format": input_format},
            "parameters": {"model": model, "tree_method": "upgma", "heterozygotes": heterozygotes,
                           "max_site_missing": max_site_missing, "min_comparable_sites": min_comparable_sites,
                           "allow_snp_correction": allow_snp_correction},
            "artifacts_sha256": artifacts, "warnings": warnings,
        })
        # mkdir is an exclusive claim: never replace an existing run, even in a race.
        output_path.mkdir()
        try:
            for p in stage.iterdir():
                p.rename(output_path / p.name)
        except BaseException:
            shutil.rmtree(output_path)
            raise
    finally:
        shutil.rmtree(stage)
    return output_path
