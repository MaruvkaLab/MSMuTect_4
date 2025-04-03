import unittest

from src.GenomicUtils.NoiseTable import get_noise_table
from src.IndelCalling import CallAllelesFast
from tests.testing_utils.generate_histograms import get_allele_histograms, real_case_histograms
from src.IndelCalling.CallMutations import *


class TestCallMutations(unittest.TestCase):

    def test_hist2vec(self):
        histograms = get_allele_histograms()
        vecs = hist2vecs(histograms[0], histograms[1])
        self.assertTrue(0 in vecs.first_set and 0 in vecs.second_set)

    def test_real(self):
        noise_table = get_noise_table()
        fisher = Fisher()
        normal, tumor = real_case_histograms()
        tumor_alleles = CallAllelesFast.calculate_alleles(tumor, noise_table,
                                                          required_read_support=5)
        normal_alleles = CallAllelesFast.calculate_alleles(normal, noise_table,
                                                           required_read_support=5)
        print(call_mutations(normal_alleles, tumor_alleles, noise_table, fisher))


if __name__ == '__main__':
    unittest.main()
