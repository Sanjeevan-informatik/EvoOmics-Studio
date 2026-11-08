"""Shared, validated domain objects, independent of input format and UI."""
from dataclasses import dataclass
import re

DNA = frozenset("ACGT")
IUPAC = frozenset("ACGTRYSWKMBDHVN-?")


class ValidationError(ValueError):
    """An input or scientific assumption prevents a valid analysis."""


@dataclass(frozen=True)
class Alignment:
    samples: tuple[str, ...]
    sequences: tuple[str, ...]
    source_kind: str = "fasta"

    def __post_init__(self):
        if len(self.samples) < 2 or len(self.samples) != len(self.sequences):
            raise ValidationError("At least two samples and one sequence per sample are required")
        if len(set(self.samples)) != len(self.samples):
            raise ValidationError("Duplicate sample IDs")
        if any(not re.fullmatch(r"[A-Za-z0-9_.-]+", s) for s in self.samples):
            raise ValidationError("Sample IDs must use only letters, digits, _, . and -")
        lengths = {len(s) for s in self.sequences}
        if len(lengths) != 1 or 0 in lengths:
            raise ValidationError("Sequences must be aligned, non-empty and equal in length")
        if any(set(s) - IUPAC for s in self.sequences):
            raise ValidationError("Unsupported DNA symbols (use uppercase DNA/IUPAC)")

    @property
    def sites(self) -> int:
        return len(self.sequences[0])
