from dataclasses import dataclass
import os
from typing import List

from src.IndelCalling.CallMutations import calculate_AICs


@dataclass
class PerMotifLengthRes:
    num_loci: List[int]
    num_mutations: List[int]

    @property
    def mutation_rate(self) -> List[int]:
        raise NotImplementedError

def get_loci_counts_per_motif_size(fp: str, out_dir: str, old: bool) -> PerMotifLengthRes:
    loci_files = [open(os.path.join(out_dir, f"{i}.txt"), 'w+') for i in range(3, 16)]
    ret = PerMotifLengthRes([0 for i in range(15)], [0 for i in range(15)])
    if old:
        call_idx = 48
    else:
        call_idx = 50
    with open(fp, 'r') as results_file:
        header_line = results_file.readline()
        current_line = results_file.readline()
        # while len(current_line.strip())!=0:
        for i in range(3_500_000):
            split_line = current_line.split("\t")
            motif = split_line[3].strip()
            call = split_line[call_idx].strip()
            motif_length = len(motif)
            ret.num_loci[motif_length-1]+=1
            if call == "M":
                ret.num_mutations[motif_length-1]+=1
                if motif_length>=3:
                    loci_files[motif_length-3].write(current_line+"\n")
            current_line = results_file.readline()

    return ret

def main():
    res_hist = get_loci_counts_per_motif_size("/home/avraham/MaruvkaLab/Texas/strict_msmutect/results/res_low_purity.full.mut.tsv",
                                              "/home/avraham/MaruvkaLab/Texas/strict_msmutect/results/new_results3", old=False)
    # res_hist = get_loci_counts_per_motif_size(
    #     "/home/avraham/MaruvkaLab/Texas/strict_msmutect/results/TCGA-BR-4201.full.mut.tsv",
    #     "/home/avraham/MaruvkaLab/Texas/strict_msmutect/results/tmp", old=True)
    print(res_hist)

if __name__ == '__main__':
    main()

