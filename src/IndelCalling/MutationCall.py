from typing import Tuple

from scipy.stats import ks_2samp
from src.IndelCalling.AlleleSet import AlleleSet
from src.IndelCalling.AICs import AICs
from src.IndelCalling.hist2vecs import hist2samps


class MutationCall:
    # pseudo enum. Not using real enums for backwards compatibility with obsolete python versions (biologists hate to update)
    LOSS_OF_HETEROZYGOSITY = -6
    REVERTED_TO_REFERENCE = -5
    NO_ALLELES = -4
    BORDERLINE_NONMUTATION = -3
    TOO_MANY_ALLELES = -2
    INSUFFICIENT = -1
    NOT_MUTATION = 0
    MUTATION = 1


    def __init__(self, call: int, normal_alleles: AlleleSet, tumor_alleles: AlleleSet, aic_values: AICs, p_value=-1):
        self.call = call
        self.normal_alleles = normal_alleles
        self.tumor_alleles = tumor_alleles
        self.aic_values = aic_values
        self.p_value = p_value

    def format_pval(self):
        if self.p_value == -1:
            return 'NA'
        else:
            return str(self.p_value)

    def call_abbreviation(self, call: int) -> str:
        abbreviations = {
                         MutationCall.LOSS_OF_HETEROZYGOSITY: "LOH",
                         MutationCall.REVERTED_TO_REFERENCE: "RR",
                         MutationCall.NO_ALLELES: "AN", # either tumor or normal lacks alleles
                         MutationCall.BORDERLINE_NONMUTATION: "FFT",  # failed fisher test
                         MutationCall.TOO_MANY_ALLELES : "TMA",
                         MutationCall.INSUFFICIENT: "INS",
                         MutationCall.NOT_MUTATION: "NM",
                         MutationCall.MUTATION: "M",
                        }
        return abbreviations[call]

    def ks_test_value(self) -> Tuple[float, float]:
        reads_samps = hist2samps(self.tumor_alleles.histogram, self.normal_alleles.histogram)  # order is important for Fisher test
        if len(reads_samps[0]) == 0 or len(reads_samps[1])== 0:
            return -1, -1
        else:
            ks_test_result = ks_2samp(reads_samps[0], reads_samps[1])
            return ks_test_result.pvalue, ks_test_result.statistic


    @staticmethod
    def header():
        return f"CALL\tFISHER_TEST_PVALUE\t{AICs.header()}\tKS_TEST_PVALUE\tKS_TEST_STATISTIC"

    def __str__(self):
        p_val, statistic = self.ks_test_value()
        return f"{self.call_abbreviation(self.call)}\t{self.format_pval()}\t{str(self.aic_values)}\t{p_val}\t{statistic}"
