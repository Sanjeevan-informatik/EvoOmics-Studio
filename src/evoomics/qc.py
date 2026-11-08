"""Callability QC. Ambiguity is never silently treated as a fifth nucleotide."""
from evoomics.models import Alignment, DNA, ValidationError


def filter_alignment(alignment: Alignment, max_site_missing: float = 1.0):
    if not 0 <= max_site_missing <= 1:
        raise ValidationError("max_site_missing must lie between 0 and 1")
    n = len(alignment.samples)
    keep = [i for i in range(alignment.sites)
            if sum(s[i] not in DNA for s in alignment.sequences) / n <= max_site_missing]
    if not keep:
        raise ValidationError("No sites remain after missingness filtering")
    filtered = Alignment(alignment.samples,
                         tuple("".join(s[i] for i in keep) for s in alignment.sequences),
                         alignment.source_kind)
    qc = {
        "input_sites": alignment.sites, "retained_sites": len(keep),
        "removed_sites": alignment.sites - len(keep),
        "samples": {name: {"callable_sites": sum(c in DNA for c in seq),
                           "missing_fraction": sum(c not in DNA for c in seq) / len(seq)}
                    for name, seq in zip(filtered.samples, filtered.sequences)},
    }
    return filtered, qc, keep
