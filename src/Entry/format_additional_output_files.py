#!/usr/bin/env python3
"""Filter MSMuTect-style .mut.tsv files by the CALL column.

Produces two outputs from one input:
  1. <stem>.<kind>.<call>.mut.tsv  - header + only rows whose CALL == call
  2. <stem>.<kind>.call_counts.tsv - tally of every CALL value seen
                                     (cross-platform `cut -f CALL | sort | uniq -c`)

Streams line by line, so memory use is flat regardless of input size.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from collections import Counter
from pathlib import Path
from typing import NamedTuple

VALID_SUFFIXES = (".full.mut.tsv", ".partial.mut.tsv")

#: MSMuTect CALL codes -> human-readable labels for the tally file.
#: Codes not listed here are passed through verbatim.
CALL_LABELS = {
    "AN": "No Alleles (AN)",
    "FFT": "Failed Fisher Test (FFT)",
    "INS": "Insufficient Support (INS)",
    "LOH": "Loss of Heterozygosity (LOH)",
    "M": "M (Mutation)",
    "NM": "Not Mutation (NM)",
    "RR": "Reversion to Reference (RR)",
    "TMA": "Too Many Alleles (TMA)"
}


class FilterResult(NamedTuple):
    """Summary of one filtering run."""

    rows_file: Path
    counts_file: Path
    kept: int
    total: int
    counts: Counter

    @property
    def dropped(self) -> int:
        return self.total - self.kept


def _open_text(path: Path):
    """Open plain or gzipped text transparently."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return open(path, "rt", newline="")


def _split_suffix(name: str) -> tuple[str, str]:
    """Return (stem, kind) for a recognised filename, else raise.

    'sampleA.full.mut.tsv' -> ('sampleA', 'full')
    """
    base = name[:-3] if name.endswith(".gz") else name
    for suffix in VALID_SUFFIXES:
        if base.endswith(suffix) and len(base) > len(suffix):
            return base[: -len(suffix)], suffix.split(".")[1]
    raise ValueError(
        f"{name!r} does not end in one of {' or '.join(VALID_SUFFIXES)} "
        "(optionally .gz). Refusing to guess the file layout."
    )


def filter_mut_calls(
    path: str | Path,
    call: str = "M",
    outdir: str | Path | None = None,
    call_column: str = "CALL",
    sort_by: str = "call",
    label_calls: bool = True,
) -> FilterResult:
    """Split a .mut.tsv file into matching rows plus a CALL tally.

    Args:
        path: input .full.mut.tsv or .partial.mut.tsv (may be .gz).
        call: CALL value to keep. Compared after stripping whitespace.
        outdir: where to write outputs. Defaults to the input's directory.
        call_column: header name of the call column.
        sort_by: 'call' to sort the tally lexicographically by code (matches
            `sort | uniq -c`), or 'count' for descending frequency.
        label_calls: expand codes to CALL_LABELS descriptions in the tally.
            Sorting always uses the raw code either way.

    Returns:
        FilterResult with output paths and the tally.

    Raises:
        ValueError: unrecognised suffix, missing column, or empty file.
    """
    path = Path(path)
    if sort_by not in {"call", "count"}:
        raise ValueError("sort_by must be 'call' or 'count'")

    stem, kind = _split_suffix(path.name)
    outdir = Path(outdir) if outdir is not None else path.parent
    outdir.mkdir(parents=True, exist_ok=True)

    # Sanitise the call value so it can't inject a path separator.
    # tag = "".join(c if c.isalnum() or c in "-_" else "_" for c in call) or "match"
    tag = "mutated_only"  # "".join(c if c.isalnum() or c in "-_" else "_" for c in call) or "match"
    rows_file = outdir / f"{stem}.{kind}.{tag}.mut.tsv"

    # rows_file = outdir / f"{stem}.{kind}.{tag}.mut.tsv"
    counts_file = outdir / f"{stem}.{kind}.call_counts.tsv"

    if rows_file.resolve() == path.resolve():
        raise ValueError(f"output would overwrite the input: {path}")

    counts: Counter = Counter()
    kept = total = 0

    with _open_text(path) as fh:
        try:
            header = next(fh).rstrip("\r\n")
        except StopIteration:
            raise ValueError(f"{path} is empty") from None

        fields = header.split("\t")
        try:
            idx = fields.index(call_column)
        except ValueError:
            raise ValueError(
                f"no {call_column!r} column in {path.name}. Found: {fields}"
            ) from None

        with open(rows_file, "w", newline="\n") as out:
            out.write(header + "\n")
            for line in fh:
                line = line.rstrip("\r\n")
                if not line:
                    continue
                row = line.split("\t")
                if len(row) <= idx:
                    raise ValueError(
                        f"{path.name} line {total + 2}: expected >{idx} columns, "
                        f"got {len(row)}"
                    )
                total += 1
                value = row[idx].strip()
                counts[value] += 1
                if value == call:
                    out.write(line + "\n")
                    kept += 1

    if sort_by == "call":
        # Sort on the raw code, not the label, so AN/FFT/INS/... stays ordered.
        ordered = sorted(counts.items())
    else:
        ordered = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))

    with open(counts_file, "w", newline="\n") as out:
        out.write("count\t" + call_column + "\n")
        for value, n in ordered:
            label = CALL_LABELS.get(value, value) if label_calls else value
            out.write(f"{n}\t{label}\n")

    return FilterResult(rows_file, counts_file, kept, total, counts)


if __name__ == '__main__':
    f="/home/avraham/MaruvkaLab/msmutect_postprocessing/data/full_gib_files/msmutect_normal0_tumor0.full.mut.tsv"
    filter_mut_calls(
        f,
        call="M",
        outdir="/home/avraham/MaruvkaLab/msmutect_development/features/output_files/",
        call_column="CALL",
        sort_by="call"
    )