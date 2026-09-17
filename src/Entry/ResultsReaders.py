from dataclasses import dataclass
from typing import List, Optional


def to_int(value: str) -> int:
    # tolerant int parse (handles "16" and "16.0")
    return int(float(value))


def str_to_bool(value: str) -> bool:
    # columns written from a python bool, e.g. KNOWN_GERMLINE_VARIANT -> "True"/"False"
    return value.strip() == "True"


def int_to_bool(value: str) -> bool:
    # legacy columns written as 0/1, e.g. the old "Noisy Locus"
    return bool(int(value))


class NamedColumns:
    """Maps a TSV header line to column indices so rows are read by name, not position.

    This keeps the readers robust to columns being added, removed, or reordered
    (e.g. adding KNOWN_GERMLINE_VARIANT or dropping "Noisy Locus"): lookups are by
    header name, and missing columns fall back to a default instead of shifting
    every downstream field.
    """

    def __init__(self, header_line: str):
        self.index = {name.strip(): i for i, name in enumerate(header_line.rstrip("\n").split("\t"))}

    def has(self, name: str) -> bool:
        return name in self.index

    def value(self, fields: List[str], name: str, dtype=str, default=None):
        if name not in self.index:
            return default
        return dtype(fields[self.index[name]].strip())

    def resolve_base(self, suffix_base: str) -> Optional[str]:
        # Return the full column base (with whatever prefix) whose '<base>_1' exists,
        # e.g. "MOTIF_REPEATS" -> "NORMAL_MOTIF_REPEATS" (or "MOTIF_REPEATS" if unprefixed).
        # Lets histogram files be read whether their columns are prefixed or not.
        tail = f"{suffix_base}_1"
        for name in self.index:
            if name == tail or name.endswith("_" + tail):
                return name[: -len("_1")]
        return None

    def group(self, fields: List[str], base: str, dtype, count: int = 6, break_str: str = "NA") -> list:
        # collect base_1, base_2, ... up to `count`, stopping at the first missing
        # column or the first break_str ("NA") value.
        ret = []
        for k in range(1, count + 1):
            name = f"{base}_{k}"
            if name not in self.index:
                break
            raw = fields[self.index[name]].strip()
            if raw == break_str:
                break
            ret.append(dtype(raw))
        return ret


@dataclass
class LocusMutationCall:
    chromosome: str
    start: int
    end: int
    pattern: str
    ref_seq: str
    num_ref_repeats: float
    normal_motif_repeats: List[float]
    normal_motif_repeat_support: List[int]
    tumor_motif_repeats: List[float]
    tumor_motif_repeat_support: List[int]
    normal_alleles: List[float]
    tumor_alleles: List[float]
    mutation_call: str
    known_germline_variant: bool = False


class ResultsReaderMutationFile:
    def __init__(self, results_file: str):
        self.results_file = open(results_file, 'r')
        self.columns = NamedColumns(next(self.results_file))

    def __iter__(self):
        return self

    def __next__(self):
        next_line = self.results_file.readline()
        if next_line == "":
            raise StopIteration
        fields = next_line.rstrip("\n").split('\t')
        if len(fields) == 1:
            raise StopIteration
        c = self.columns
        return LocusMutationCall(
            chromosome=c.value(fields, "CHROMOSOME", str),
            start=c.value(fields, "START", int),
            end=c.value(fields, "END", int),
            pattern=c.value(fields, "PATTERN", str),
            ref_seq=c.value(fields, "REFERENCE_SEQUENCE", str),
            num_ref_repeats=c.value(fields, "REFERENCE_REPEATS", float),
            normal_motif_repeats=c.group(fields, "NORMAL_MOTIF_REPEATS", float),
            normal_motif_repeat_support=c.group(fields, "NORMAL_SUPPORTING_READS", to_int),
            tumor_motif_repeats=c.group(fields, "TUMOR_MOTIF_REPEATS", float),
            tumor_motif_repeat_support=c.group(fields, "TUMOR_SUPPORTING_READS", to_int),
            normal_alleles=c.group(fields, "NORMAL_ALLELE", float, count=4),
            tumor_alleles=c.group(fields, "TUMOR_ALLELE", float, count=4),
            mutation_call=c.value(fields, "CALL", str),
            known_germline_variant=c.value(fields, "KNOWN_GERMLINE_VARIANT", str_to_bool, default=False),
        )

    def __del__(self):
        try:
            self.results_file.close()
        except AttributeError:
            pass
