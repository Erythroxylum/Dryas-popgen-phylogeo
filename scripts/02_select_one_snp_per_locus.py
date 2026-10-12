#!/usr/bin/env python3
"""Select the best-called SNP from each ipyrad locus in a VCF.

The script writes a two-column CHROM/POS regions file for bcftools and a table
documenting the selected record. Ties are resolved by input order.
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("rt")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", type=Path)
    parser.add_argument("regions", type=Path)
    parser.add_argument("report", type=Path)
    args = parser.parse_args()

    best: dict[str, tuple[int, int, str, str, int, int]] = {}
    order = 0
    with open_text(args.vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                raise SystemExit("VCF record has no sample columns")
            chrom, pos, marker = fields[0], int(fields[1]), fields[2]
            match = re.match(r"^(loc\d+)_pos\d+$", marker)
            if not match:
                raise SystemExit(f"Marker ID is not locN_posN: {marker}")
            locus = match.group(1)
            fmt = fields[8].split(":")
            if "GT" not in fmt:
                raise SystemExit(f"GT is absent at {chrom}:{pos}")
            gt_index = fmt.index("GT")
            called = 0
            for sample in fields[9:]:
                values = sample.split(":")
                gt = values[gt_index] if gt_index < len(values) else "."
                if gt not in {".", "./.", ".|."} and "." not in gt:
                    called += 1
            candidate = (called, -order, chrom, marker, pos, order)
            current = best.get(locus)
            if current is None or candidate[:2] > current[:2]:
                best[locus] = candidate
            order += 1

    selected = sorted(best.items(), key=lambda item: item[1][5])
    args.regions.parent.mkdir(parents=True, exist_ok=True)
    with args.regions.open("w") as out:
        for _, (_, _, chrom, _, pos, _) in selected:
            out.write(f"{chrom}\t{pos}\n")
    with args.report.open("w") as out:
        out.write("locus\tCHROM\tPOS\tID\tcalled_samples\n")
        for locus, (called, _, chrom, marker, pos, _) in selected:
            out.write(f"{locus}\t{chrom}\t{pos}\t{marker}\t{called}\n")
    print(f"Selected {len(selected)} SNPs from {len(best)} ipyrad loci")


if __name__ == "__main__":
    main()

