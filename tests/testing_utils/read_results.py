from typing import List, Optional
from dataclasses import dataclass

from src.Entry.ResultsReaders import (
    NamedColumns, to_int, str_to_bool, int_to_bool,
    LocusMutationCall, ResultsReaderMutationFile,
)


@dataclass
class ResultsLine:
    chromosome: str
    start: int
    end: int
    pattern: str
    ref_seq: str
    num_ref_repeats: float
    motif_repeats: List[float]
    motif_repeat_support: List[int]
    noisy: Optional[bool] = None            # "Noisy Locus" column (absent in current output)
    known_germline_variant: bool = False


class ResultsReader:
    """Reads a single-sample histogram (.hist.tsv) file by column name."""

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
        return ResultsLine(
            chromosome=c.value(fields, "CHROMOSOME", str),
            start=c.value(fields, "START", int),
            end=c.value(fields, "END", int),
            pattern=c.value(fields, "PATTERN", str),
            ref_seq=c.value(fields, "REFERENCE_SEQUENCE", str),
            num_ref_repeats=c.value(fields, "REFERENCE_REPEATS", float),
            motif_repeats=c.group(fields, "MOTIF_REPEATS", float),
            motif_repeat_support=c.group(fields, "SUPPORTING_READS", to_int),
            noisy=c.value(fields, "Noisy Locus", int_to_bool, default=None),
            known_germline_variant=c.value(fields, "KNOWN_GERMLINE_VARIANT", str_to_bool, default=False),
        )

    def __del__(self):
        try:
            self.results_file.close()
        except AttributeError:
            pass


if __name__ == '__main__':
    a = ResultsReader("/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/test_results/mapping.hist.tsv")
    for i in range(5):
        cur_line = next(a)
        print(cur_line.motif_repeats)
        print(cur_line.motif_repeat_support)
