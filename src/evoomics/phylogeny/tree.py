"""Deterministic UPGMA with cluster-size weighted average distances."""
import math
from evoomics.models import ValidationError


def upgma(samples, matrix):
    n = len(samples)
    if n < 2 or len(set(samples)) != n or len(matrix) != n or any(len(row) != n for row in matrix):
        raise ValidationError("Invalid distance matrix dimensions or sample IDs")
    for i in range(n):
        for j in range(n):
            d = matrix[i][j]
            if not math.isfinite(d) or d < 0 or not math.isclose(d, matrix[j][i], abs_tol=1e-12):
                raise ValidationError("Distances must be finite, nonnegative and symmetric")
        if matrix[i][i] != 0:
            raise ValidationError("Distance diagonal must be zero")
    # Quote Newick labels, including escaping embedded single quotes.
    clusters = {i: ("'" + s.replace("'", "''") + "'", 1, 0.0) for i, s in enumerate(samples)}
    distances = {(i, j): matrix[i][j] for i in range(n) for j in range(i + 1, n)}
    next_id = n
    while len(clusters) > 1:
        a, b = min(distances, key=lambda pair: (distances[pair], pair))
        height = distances[(a, b)] / 2
        tree_a, size_a, height_a = clusters[a]
        tree_b, size_b, height_b = clusters[b]
        if height + 1e-12 < max(height_a, height_b):
            raise ValidationError("Non-monotonic UPGMA merge")
        tree = f"({tree_a}:{max(0, height-height_a):.10g},{tree_b}:{max(0, height-height_b):.10g})"
        updates = {}
        for k in clusters:
            if k not in (a, b):
                da = distances[tuple(sorted((a, k)))]
                db = distances[tuple(sorted((b, k)))]
                updates[(k, next_id)] = (size_a * da + size_b * db) / (size_a + size_b)
        distances = {pair: d for pair, d in distances.items() if a not in pair and b not in pair}
        distances.update(updates)
        del clusters[a], clusters[b]
        clusters[next_id] = tree, size_a + size_b, height
        next_id += 1
    return next(iter(clusters.values()))[0] + ";"
