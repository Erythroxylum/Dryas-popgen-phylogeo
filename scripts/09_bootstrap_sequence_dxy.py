#!/usr/bin/env python3
"""Bootstrap filtered sequence dXY by resampling complete ipyrad loci."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

import numpy as np


def write_matrix(path: Path, populations: list[str], pair_i, pair_j, values) -> None:
    matrix = np.full((len(populations), len(populations)), np.nan)
    matrix[pair_i, pair_j] = values
    matrix[pair_j, pair_i] = values
    np.fill_diagonal(matrix, 0)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["population", *populations])
        for i, population in enumerate(populations):
            writer.writerow([
                population,
                *("" if not np.isfinite(x) else f"{x:.12g}" for x in matrix[i]),
            ])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("per_locus_npz", type=Path)
    parser.add_argument("output_prefix", type=Path)
    parser.add_argument("--nboot", type=int, default=int(os.getenv("NBOOT", "1000")))
    parser.add_argument("--seed", type=int, default=int(os.getenv("BOOTSTRAP_SEED", "20261012")))
    parser.add_argument("--batch-size", type=int, default=25)
    args = parser.parse_args()
    if args.nboot < 2:
        raise SystemExit("Use at least two bootstrap replicates")

    dat = np.load(args.per_locus_npz, allow_pickle=False)
    differences = dat["differences"]
    comparisons = dat["comparisons"]
    pair_i, pair_j = dat["pair_i"], dat["pair_j"]
    populations = dat["populations"].astype(str).tolist()
    classes = dat["pop_class"].astype(str).tolist()
    n_loci, n_pairs = differences.shape
    rng = np.random.default_rng(args.seed)
    estimates = np.empty((args.nboot, n_pairs), dtype=np.float64)

    print(
        f"Bootstrapping {n_pairs} population pairs across {n_loci} loci "
        f"for {args.nboot} replicates...", flush=True
    )
    diff_float = differences.astype(np.float64)
    comp_float = comparisons.astype(np.float64)
    for start in range(0, args.nboot, args.batch_size):
        stop = min(start + args.batch_size, args.nboot)
        weights = np.zeros((stop - start, n_loci), dtype=np.float64)
        for row in range(stop - start):
            sampled = rng.integers(0, n_loci, size=n_loci)
            weights[row] = np.bincount(sampled, minlength=n_loci)
        numerator = weights @ diff_float
        denominator = weights @ comp_float
        np.divide(
            numerator, denominator, out=estimates[start:stop],
            where=denominator > 0
        )
        estimates[start:stop][denominator <= 0] = np.nan
        print(f"Completed {stop}/{args.nboot} bootstraps", flush=True)

    total_diff = differences.sum(axis=0, dtype=np.uint64)
    total_comp = comparisons.sum(axis=0, dtype=np.uint64)
    point = np.divide(
        total_diff, total_comp,
        out=np.full(n_pairs, np.nan), where=total_comp > 0
    )
    lower = np.nanquantile(estimates, 0.025, axis=0)
    upper = np.nanquantile(estimates, 0.975, axis=0)
    n_loci_pair = (comparisons > 0).sum(axis=0)

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    long_path = Path(str(args.output_prefix) + "_bootstrap_ci.tsv")
    with long_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "population_1", "population_2", "class_1", "class_2",
            "n_loci", "dxy", "lower_95", "upper_95", "nboot"
        ])
        for k, (i, j) in enumerate(zip(pair_i, pair_j)):
            writer.writerow([
                populations[i], populations[j], classes[i], classes[j],
                int(n_loci_pair[k]), point[k], lower[k], upper[k], args.nboot
            ])

    write_matrix(Path(str(args.output_prefix) + "_matrix.tsv"),
                 populations, pair_i, pair_j, point)
    write_matrix(Path(str(args.output_prefix) + "_lower95_matrix.tsv"),
                 populations, pair_i, pair_j, lower)
    write_matrix(Path(str(args.output_prefix) + "_upper95_matrix.tsv"),
                 populations, pair_i, pair_j, upper)
    print(f"Wrote {long_path}")


if __name__ == "__main__":
    main()
