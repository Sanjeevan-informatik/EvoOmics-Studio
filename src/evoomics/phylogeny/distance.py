"""Pairwise deletion with reported comparable-site counts."""
import math
from evoomics.models import Alignment, DNA, ValidationError


MODELS = ("p", "jc69", "k2p")
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def pair_distance(a: str, b: str, model: str = "p", min_sites: int = 1):
    if model not in MODELS or min_sites < 1 or len(a) != len(b):
        raise ValidationError("Invalid distance model, minimum sites or sequence lengths")
    count = differences = transitions = 0
    for x, y in zip(a, b):
        if x not in DNA or y not in DNA:
            continue
        count += 1
        differences += x != y
        transitions += (x, y) in TRANSITIONS
    if count < min_sites:
        raise ValidationError(f"Only {count} comparable sites; need {min_sites}")
    p = differences / count
    if model == "p":
        return p, count
    if model == "jc69":
        term = 1 - 4 * p / 3
        if term <= 0:
            raise ValidationError("JC69 saturated: mismatch fraction >= 0.75")
        return -0.75 * math.log(term), count
    transition_fraction = transitions / count
    transversion_fraction = (differences - transitions) / count
    first, second = 1 - 2 * transition_fraction - transversion_fraction, 1 - 2 * transversion_fraction
    if first <= 0 or second <= 0:
        raise ValidationError("K2P saturated: logarithm outside its domain")
    return -0.5 * math.log(first) - 0.25 * math.log(second), count


def distance_matrix(alignment: Alignment, model: str = "p", min_sites: int = 1):
    n = len(alignment.samples)
    matrix = [[0.0] * n for _ in range(n)]
    counts = [[0] * n for _ in range(n)]
    for i, seq in enumerate(alignment.sequences):
        counts[i][i] = sum(c in DNA for c in seq)
        for j in range(i):
            try:
                distance, count = pair_distance(seq, alignment.sequences[j], model, min_sites)
            except ValidationError as exc:
                raise ValidationError(f"{alignment.samples[i]} vs {alignment.samples[j]}: {exc}") from exc
            matrix[i][j] = matrix[j][i] = distance
            counts[i][j] = counts[j][i] = count
    return matrix, counts
