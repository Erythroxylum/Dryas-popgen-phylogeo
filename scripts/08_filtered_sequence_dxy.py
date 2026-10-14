#!/usr/bin/env python3
"""Calculate callable-site-weighted pairwise dXY from filtered ipyrad loci.

Invariant callable positions come from the ipyrad .loci alignment. Differences
are introduced only at biallelic variants retained in the cleaned VCF, so
masked genotypes and batch-filtered loci cannot re-enter through raw consensus
sequences. A compressed locus-by-population-pair file is retained for locus
bootstrapping by script 09.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path

import numpy as np


VALID = np.zeros(256, dtype=np.uint8)
for value in b"ACGTRYSWKMBDHV":
    VALID[value] = 1


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("rt")


def load_aliases(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            previous = result.get(row["raw_id"])
            if previous is not None and previous != row["sampleID"]:
                raise SystemExit(f"Ambiguous alias: {row['raw_id']}")
            result[row["raw_id"]] = row["sampleID"]
    return result


def load_popmap(path: Path):
    sample_to_pop: dict[str, str] = {}
    pop_class: dict[str, str] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["selected"].lower() not in {"true", "1"}:
                continue
            population = row["pop_gen_population"]
            sample_to_pop[row["sampleID"]] = population
            old = pop_class.get(population)
            if old is not None and old != row["pop_gen"]:
                raise SystemExit(f"Population has mixed classes: {population}")
            pop_class[population] = row["pop_gen"]
    populations = sorted(pop_class)
    return sample_to_pop, populations, pop_class


def load_excluded(path: Path | None) -> set[int]:
    if path is None:
        return set()
    excluded: set[int] = set()
    with path.open() as handle:
        for line in handle:
            token = line.strip().split()[0] if line.strip() else ""
            if not token:
                continue
            match = re.search(r"(?:loc)?(\d+)", token)
            if not match:
                raise SystemExit(f"Cannot parse excluded locus: {token}")
            excluded.add(int(match.group(1)))
    return excluded


def parse_gt(value: str) -> int | None:
    if value in {".", "./.", ".|."} or "." in value:
        return None
    alleles = re.split(r"[/|]", value)
    if len(alleles) != 2 or any(x not in {"0", "1"} for x in alleles):
        return None
    return sum(int(x) for x in alleles)


def load_vcf_summaries(vcf, aliases, sample_to_pop, populations, excluded):
    pop_index = {p: i for i, p in enumerate(populations)}
    by_locus: dict[int, dict[int, tuple[np.ndarray, np.ndarray]]] = defaultdict(dict)
    sample_pops: list[int | None] = []
    n_sites = 0
    with open_text(vcf) as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                raw_ids = line.rstrip("\n").split("\t")[9:]
                canonical = [aliases.get(x, x) for x in raw_ids]
                sample_pops = [
                    pop_index[sample_to_pop[x]] if x in sample_to_pop else None
                    for x in canonical
                ]
                continue
            if line.startswith("#"):
                continue
            if not sample_pops:
                raise SystemExit("VCF #CHROM header not found")
            fields = line.rstrip("\n").split("\t")
            match = re.match(r"^loc(\d+)_pos(\d+)$", fields[2])
            if not match:
                raise SystemExit(f"Marker ID is not locN_posN: {fields[2]}")
            locus, position = int(match.group(1)), int(match.group(2))
            if locus in excluded:
                continue
            fmt = fields[8].split(":")
            if "GT" not in fmt:
                raise SystemExit(f"GT absent for {fields[2]}")
            gt_index = fmt.index("GT")
            alt = np.zeros(len(populations), dtype=np.int64)
            chromosomes = np.zeros(len(populations), dtype=np.int64)
            for p, sample in zip(sample_pops, fields[9:]):
                if p is None:
                    continue
                values = sample.split(":")
                gt = values[gt_index] if gt_index < len(values) else "."
                dosage = parse_gt(gt)
                if dosage is not None:
                    alt[p] += dosage
                    chromosomes[p] += 2
            by_locus[locus][position] = (alt, chromosomes)
            n_sites += 1
    return by_locus, n_sites


def write_matrix(path: Path, populations: list[str], values: np.ndarray) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["population", *populations])
        for i, population in enumerate(populations):
            writer.writerow([
                population,
                *("" if not np.isfinite(x) else f"{x:.12g}" for x in values[i]),
            ])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("loci", type=Path)
    parser.add_argument("clean_vcf", type=Path)
    parser.add_argument("sample_map", type=Path)
    parser.add_argument("aliases", type=Path)
    parser.add_argument("output_prefix", type=Path)
    parser.add_argument("--exclude-loci", type=Path)
    parser.add_argument("--position-base", type=int, choices=(0, 1), default=0)
    args = parser.parse_args()

    aliases = load_aliases(args.aliases)
    sample_to_pop, populations, pop_class = load_popmap(args.sample_map)
    excluded = load_excluded(args.exclude_loci)
    pop_index = {p: i for i, p in enumerate(populations)}
    pair_i, pair_j = np.triu_indices(len(populations), 1)

    print("Reading cleaned VCF genotype summaries...", flush=True)
    retained, n_sites = load_vcf_summaries(
        args.clean_vcf, aliases, sample_to_pop, populations, excluded
    )
    print(f"Loaded {n_sites} retained variant sites", flush=True)

    locus_ids: list[int] = []
    locus_differences: list[np.ndarray] = []
    locus_comparisons: list[np.ndarray] = []
    samples_seen: set[str] = set()

    def process_locus(locus: int, sequences: dict[str, str]) -> None:
        if locus in excluded or not sequences:
            return
        canonical: dict[str, str] = {}
        for raw_id, sequence in sequences.items():
            sample = aliases.get(raw_id, raw_id)
            if sample in sample_to_pop:
                canonical[sample] = sequence.upper()
                samples_seen.add(sample)
        if not canonical:
            return
        length = max(map(len, canonical.values()))
        called = np.zeros((length, len(populations)), dtype=np.int16)
        for sample, sequence in canonical.items():
            encoded = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
            called[: len(encoded), pop_index[sample_to_pop[sample]]] += 2 * VALID[encoded]

        called64 = called.astype(np.int64, copy=False)
        comparisons = called64.T @ called64
        differences = np.zeros_like(comparisons)
        for position, (alt, chromosomes) in retained.get(locus, {}).items():
            zero_position = position - args.position_base
            raw = (called64[zero_position] if 0 <= zero_position < length
                   else np.zeros(len(populations), dtype=np.int64))
            comparisons += np.outer(chromosomes, chromosomes) - np.outer(raw, raw)
            ref = chromosomes - alt
            differences += np.outer(alt, ref) + np.outer(ref, alt)

        pair_comp = comparisons[pair_i, pair_j]
        if np.any(pair_comp > 0):
            if np.any(pair_comp < 0):
                raise SystemExit(f"Negative callable comparison count at locus {locus}")
            locus_ids.append(locus)
            locus_comparisons.append(pair_comp.astype(np.uint32))
            locus_differences.append(differences[pair_i, pair_j].astype(np.uint32))

    current: dict[str, str] = {}
    locus = 0
    with open_text(args.loci) as handle:
        for line in handle:
            stripped = line.strip()
            if stripped.startswith("//"):
                process_locus(locus, current)
                locus += 1
                current = {}
                if locus % 1000 == 0:
                    print(f"Processed {locus} loci", flush=True)
                continue
            if not stripped:
                continue
            fields = stripped.split()
            if len(fields) >= 2:
                raw_id = fields[0].lstrip(">")
                if raw_id in aliases or raw_id in sample_to_pop:
                    current[raw_id] = fields[1]
    if current:
        process_locus(locus, current)

    missing = sorted(set(sample_to_pop) - samples_seen)
    if missing:
        raise SystemExit(
            "Selected samples absent from .loci after alias resolution: "
            + ", ".join(missing)
        )
    if not locus_ids:
        raise SystemExit("No usable loci were found")

    differences = np.stack(locus_differences)
    comparisons = np.stack(locus_comparisons)
    total_diff = differences.sum(axis=0, dtype=np.uint64)
    total_comp = comparisons.sum(axis=0, dtype=np.uint64)
    dxy = np.divide(
        total_diff, total_comp,
        out=np.full(len(pair_i), np.nan), where=total_comp > 0
    )
    n_loci_pair = (comparisons > 0).sum(axis=0)

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    long_path = Path(str(args.output_prefix) + ".tsv")
    matrix_path = Path(str(args.output_prefix) + "_matrix.tsv")
    npz_path = Path(str(args.output_prefix) + "_per_locus.npz")
    with long_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "population_1", "population_2", "class_1", "class_2", "n_loci",
            "different_allele_pairs", "callable_allele_pairs", "dxy"
        ])
        for k, (i, j) in enumerate(zip(pair_i, pair_j)):
            writer.writerow([
                populations[i], populations[j], pop_class[populations[i]],
                pop_class[populations[j]], int(n_loci_pair[k]),
                int(total_diff[k]), int(total_comp[k]), dxy[k]
            ])

    matrix = np.full((len(populations), len(populations)), np.nan)
    matrix[pair_i, pair_j] = dxy
    matrix[pair_j, pair_i] = dxy
    np.fill_diagonal(matrix, 0)
    write_matrix(matrix_path, populations, matrix)
    np.savez_compressed(
        npz_path, locus=np.asarray(locus_ids, dtype=np.int32),
        differences=differences, comparisons=comparisons,
        pair_i=pair_i.astype(np.int16), pair_j=pair_j.astype(np.int16),
        populations=np.asarray(populations),
        pop_class=np.asarray([pop_class[x] for x in populations])
    )
    print(f"Wrote {long_path}")
    print(f"Wrote {matrix_path}")
    print(f"Wrote {npz_path}")


if __name__ == "__main__":
    main()
