#!/usr/bin/env python3
"""Calculate sequence-length-normalized nucleotide diversity by population.

The ipyrad .loci file supplies callable sequence positions and invariant sites.
Only variants retained in the cleaned VCF contribute nucleotide differences.
Thus, masked or filtered genotypes do not re-enter through raw IUPAC consensus
sequences. The VCF may be in genomic order; records are keyed by locN_posN.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path


VALID_BASES = set("ACGTRYSWKMBDHV")


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("rt")


def load_aliases(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            old = mapping.get(row["raw_id"])
            if old is not None and old != row["sampleID"]:
                raise SystemExit(f"Ambiguous alias: {row['raw_id']}")
            mapping[row["raw_id"]] = row["sampleID"]
    return mapping


def load_popmap(path: Path) -> tuple[dict[str, str], list[str], dict[str, str]]:
    sample_to_pop: dict[str, str] = {}
    sample_class: dict[str, str] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["selected"].lower() not in {"true", "1"}:
                continue
            sample_to_pop[row["sampleID"]] = row["pop_gen_population"]
            sample_class[row["sampleID"]] = row["pop_gen"]
    populations = sorted(set(sample_to_pop.values()))
    return sample_to_pop, populations, sample_class


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
    if len(alleles) != 2 or any(a not in {"0", "1"} for a in alleles):
        return None
    return sum(int(a) for a in alleles)


def load_vcf_summaries(
    vcf: Path,
    aliases: dict[str, str],
    sample_to_pop: dict[str, str],
    populations: list[str],
    excluded: set[int],
) -> dict[tuple[int, int], tuple[bytes, bytes]]:
    pop_index = {p: i for i, p in enumerate(populations)}
    summaries: dict[tuple[int, int], tuple[bytes, bytes]] = {}
    sample_pops: list[int | None] = []
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
            marker = fields[2]
            match = re.match(r"^loc(\d+)_pos(\d+)$", marker)
            if not match:
                raise SystemExit(f"Marker ID is not locN_posN: {marker}")
            locus, position = int(match.group(1)), int(match.group(2))
            if locus in excluded:
                continue
            fmt = fields[8].split(":")
            if "GT" not in fmt:
                raise SystemExit(f"GT absent for {marker}")
            gt_index = fmt.index("GT")
            alt = bytearray(len(populations))
            chromosomes = bytearray(len(populations))
            for pop_i, sample in zip(sample_pops, fields[9:]):
                if pop_i is None:
                    continue
                values = sample.split(":")
                gt = values[gt_index] if gt_index < len(values) else "."
                dosage = parse_gt(gt)
                if dosage is None:
                    continue
                alt[pop_i] += dosage
                chromosomes[pop_i] += 2
            summaries[(locus, position)] = (bytes(alt), bytes(chromosomes))
    return summaries


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
    sample_to_pop, populations, sample_class = load_popmap(args.sample_map)
    pop_class = {}
    for sample, population in sample_to_pop.items():
        cls = sample_class[sample]
        if population in pop_class and pop_class[population] != cls:
            raise SystemExit(f"Population has mixed classes: {population}")
        pop_class[population] = cls
    excluded = load_excluded(args.exclude_loci)
    print("Reading cleaned VCF genotype summaries...", flush=True)
    site_summary = load_vcf_summaries(
        args.clean_vcf, aliases, sample_to_pop, populations, excluded
    )
    print(f"Loaded {len(site_summary)} retained variant sites", flush=True)

    pop_index = {p: i for i, p in enumerate(populations)}
    total_differences = [0] * len(populations)
    total_comparisons = [0] * len(populations)
    loci_used = [0] * len(populations)
    sites_used = [0] * len(populations)
    samples_seen: set[str] = set()
    locus_rows: list[tuple] = []

    def process_locus(locus_number: int, sequences: dict[str, str]) -> None:
        if locus_number in excluded or not sequences:
            return
        canonical_sequences: dict[str, str] = {}
        for raw_id, sequence in sequences.items():
            sample = aliases.get(raw_id, raw_id)
            if sample in sample_to_pop:
                canonical_sequences[sample] = sequence.upper()
                samples_seen.add(sample)
        if not canonical_sequences:
            return
        length = max(map(len, canonical_sequences.values()))
        locus_diff = [0] * len(populations)
        locus_comp = [0] * len(populations)
        locus_sites = [0] * len(populations)
        for zero_pos in range(length):
            vcf_pos = zero_pos + args.position_base
            retained = site_summary.get((locus_number, vcf_pos))
            if retained is not None:
                alt, chromosomes = retained
                for p in range(len(populations)):
                    n = chromosomes[p]
                    if n >= 2:
                        locus_comp[p] += n * (n - 1) // 2
                        locus_diff[p] += alt[p] * (n - alt[p])
                        locus_sites[p] += 1
                continue
            called_individuals = [0] * len(populations)
            for sample, sequence in canonical_sequences.items():
                if zero_pos < len(sequence) and sequence[zero_pos] in VALID_BASES:
                    called_individuals[pop_index[sample_to_pop[sample]]] += 1
            for p, n_individuals in enumerate(called_individuals):
                n = 2 * n_individuals
                if n >= 2:
                    locus_comp[p] += n * (n - 1) // 2
                    locus_sites[p] += 1
        for p, population in enumerate(populations):
            if locus_comp[p] > 0:
                loci_used[p] += 1
                total_differences[p] += locus_diff[p]
                total_comparisons[p] += locus_comp[p]
                sites_used[p] += locus_sites[p]
                locus_rows.append(
                    (locus_number, population, pop_class[population],
                     locus_diff[p], locus_comp[p], locus_sites[p])
                )

    current_sequences: dict[str, str] = {}
    locus_number = 0
    with open_text(args.loci) as handle:
        for line in handle:
            stripped = line.strip()
            if stripped.startswith("//"):
                process_locus(locus_number, current_sequences)
                locus_number += 1
                current_sequences = {}
                if locus_number % 1000 == 0:
                    print(f"Processed {locus_number} loci", flush=True)
                continue
            if not stripped:
                continue
            fields = stripped.split()
            if len(fields) >= 2:
                raw_id = fields[0].lstrip(">")
                sequence = fields[1]
                if raw_id in aliases or raw_id in sample_to_pop:
                    current_sequences[raw_id] = sequence
    if current_sequences:
        process_locus(locus_number, current_sequences)
        locus_number += 1

    missing_selected = sorted(set(sample_to_pop) - samples_seen)
    if missing_selected:
        raise SystemExit(
            "Selected sample IDs were never observed in .loci after alias resolution: "
            + ", ".join(missing_selected)
        )

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    population_file = Path(str(args.output_prefix) + "_population_pi.tsv")
    locus_file = Path(str(args.output_prefix) + "_per_locus_pi.tsv.gz")
    with population_file.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "population", "pop_gen", "n_loci", "callable_sites",
            "different_allele_pairs", "callable_allele_pairs", "pi"
        ])
        for p, population in enumerate(populations):
            pi = (total_differences[p] / total_comparisons[p]
                  if total_comparisons[p] else "")
            writer.writerow([
                population, pop_class[population], loci_used[p], sites_used[p],
                total_differences[p], total_comparisons[p], pi
            ])
    with gzip.open(locus_file, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "locus", "population", "pop_gen", "different_allele_pairs",
            "callable_allele_pairs", "callable_sites"
        ])
        writer.writerows(locus_rows)
    print(f"Wrote {population_file}")
    print(f"Wrote {locus_file}")


if __name__ == "__main__":
    main()

