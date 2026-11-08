"""Small-cohort FASTA and VCF readers. No external commands are executed."""
import gzip
from pathlib import Path
import re

from evoomics.models import Alignment, DNA, ValidationError


LIMIT_SAMPLES = 200
LIMIT_SITES = 100_000


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def read_fasta(path: Path) -> tuple[Alignment, dict]:
    names, seqs = [], []
    lengths = []
    with open_text(path) as handle:
        for number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if not line[1:].strip():
                    raise ValidationError(f"FASTA line {number}: empty header")
                names.append(line[1:].split()[0])
                seqs.append([])
                lengths.append(0)
                if len(names) > LIMIT_SAMPLES:
                    raise ValidationError(f"Foundation limit: {LIMIT_SAMPLES} samples")
            else:
                if not names:
                    raise ValidationError(f"FASTA line {number}: sequence before header")
                lengths[-1] += len(line)
                if lengths[-1] > LIMIT_SITES:
                    raise ValidationError(f"Foundation limit: {LIMIT_SITES} sites")
                seqs[-1].append(line.upper())
    alignment = Alignment(tuple(names), tuple("".join(s) for s in seqs))
    return alignment, {"input_sites": alignment.sites}


IUPAC_CODES = {
    frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("CG"): "S",
    frozenset("AT"): "W", frozenset("GT"): "K", frozenset("AC"): "M",
}


def read_vcf(path: Path, heterozygotes: str = "mask") -> tuple[Alignment, dict]:
    """Read PASS biallelic SNPs; accept haploid/diploid GT, mask missing calls.

    A VCF genotype is not a phased haplotype. IUPAC calls are kept in exports
    but excluded from distance denominators. Unfiltered '.' records are skipped.
    """
    if heterozygotes not in {"mask", "iupac"}:
        raise ValidationError("Heterozygote policy must be mask or iupac")
    names, sequences = None, []
    stats = dict(records=0, retained_snps=0, skipped_filter=0, skipped_non_biallelic_snp=0,
                 missing_genotypes=0, heterozygous_genotypes=0)
    seen = set()
    with open_text(path) as handle:
        for number, raw in enumerate(handle, 1):
            line = raw.rstrip("\r\n")
            if not line or line.startswith("##"):
                continue
            fields = line.split("\t")
            if line.startswith("#CHROM"):
                if names is not None:
                    raise ValidationError("Duplicate VCF header")
                expected = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
                if fields[:9] != expected or len(fields) < 11:
                    raise ValidationError("VCF needs standard columns and at least two samples")
                names = fields[9:]
                if len(names) > LIMIT_SAMPLES:
                    raise ValidationError(f"Foundation limit: {LIMIT_SAMPLES} samples")
                sequences = [[] for _ in names]
                continue
            if line.startswith("#"):
                continue
            if names is None:
                raise ValidationError("VCF record before #CHROM header")
            if len(fields) != 9 + len(names):
                raise ValidationError(f"VCF line {number}: incorrect column count")
            stats["records"] += 1
            try:
                pos = int(fields[1])
                if pos < 1 or not fields[0]:
                    raise ValueError
            except ValueError as exc:
                raise ValidationError(f"VCF line {number}: invalid coordinate") from exc
            if fields[6] != "PASS":
                stats["skipped_filter"] += 1
                continue
            ref, alt = fields[3:5]
            if ref not in DNA or alt not in DNA or ref == alt:
                stats["skipped_non_biallelic_snp"] += 1
                continue
            coordinate = (fields[0], pos)
            if coordinate in seen:
                raise ValidationError(f"Duplicate retained VCF coordinate: {coordinate}")
            seen.add(coordinate)
            formats = fields[8].split(":")
            if len(set(formats)) != len(formats) or "GT" not in formats:
                raise ValidationError(f"VCF line {number}: missing or duplicate FORMAT keys")
            gt_index = formats.index("GT")
            for i, cell in enumerate(fields[9:]):
                values = cell.split(":")
                gt = values[gt_index] if gt_index < len(values) else "."
                if not re.fullmatch(r"[01.](?:[/|][01.])?", gt):
                    raise ValidationError(f"VCF line {number}: unsupported GT {gt!r}")
                alleles = re.split(r"[/|]", gt)
                if "." in alleles:
                    base = "N"
                    stats["missing_genotypes"] += 1
                elif len(set(alleles)) > 1:
                    stats["heterozygous_genotypes"] += 1
                    base = IUPAC_CODES[frozenset((ref, alt))] if heterozygotes == "iupac" else "N"
                else:
                    base = (ref, alt)[int(alleles[0])]
                sequences[i].append(base)
            stats["retained_snps"] += 1
            if stats["retained_snps"] > LIMIT_SITES:
                raise ValidationError(f"Foundation limit: {LIMIT_SITES} retained sites")
    if names is None:
        raise ValidationError("Missing VCF #CHROM header")
    return Alignment(tuple(names), tuple("".join(s) for s in sequences), "vcf"), stats
