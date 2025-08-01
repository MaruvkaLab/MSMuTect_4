from dataclasses import dataclass
from typing import List




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


class ResultsReaderMutationFile:
    def __init__(self, results_file: str):
        self.results_file = open(results_file, 'r')
        header_line = next(self.results_file)

    def __iter__(self):
        return self

    def variable_length_tsv_list(self, broken_line: List[str], start: int, end: int, dtype: type, break_str="NA"):
        ret = []
        for i in range(start, end):
            if broken_line[i].strip() == break_str:
                break
            else:
                ret.append(dtype(broken_line[i]))
        return ret

    def __next__(self):
        next_line = self.results_file.readline()
        if next_line == "":
            raise StopIteration
        else:
            broken_line = next_line.split('\t')
            if len(broken_line) == 1:
                raise StopIteration
            dtypes = [str, int, int, str, str, float]
            chromosome, start, end, pattern, ref_seq, num_ref_repeats = ((dtypes[i](broken_line[i]) for i in range(len(dtypes))))

            normal_motif_repeats = self.variable_length_tsv_list(broken_line, 6, 12, float)
            normal_motif_repeat_support = self.variable_length_tsv_list(broken_line, 12, 18, int)

            tumor_motif_repeats = self.variable_length_tsv_list(broken_line, 27, 33, float)
            tumor_motif_repeats_support = self.variable_length_tsv_list(broken_line, 33, 39, int)

            normal_alleles = self.variable_length_tsv_list(broken_line, 19, 23, float)
            tumor_alleles = self.variable_length_tsv_list(broken_line, 41, 45, float)

            mutation_call = broken_line[48].strip()
            return LocusMutationCall(chromosome, start, end, pattern, ref_seq, num_ref_repeats, normal_motif_repeats, normal_motif_repeat_support,
                                           tumor_motif_repeats, tumor_motif_repeats_support, normal_alleles, tumor_alleles, mutation_call)

    def __del__(self):
        self.results_file.close()
