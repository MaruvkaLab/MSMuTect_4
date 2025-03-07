import unittest, numpy as np

from src.GenomicUtils.NoiseTable import get_noise_table
from tests.testing_utils.generate_histograms import get_allele_histograms
from src.IndelCalling.CallAlleles import *


class TestCalcAlleles(unittest.TestCase):

    def test_calculate_alleles_wminimal_support(self):
        histograms = get_allele_histograms()
        noise_table = get_noise_table()
        alleles_0 = calculate_alleles(histograms[0], noise_table, 5)
        self.assertEqual(len(alleles_0.repeat_lengths), 2)

    def test_calculate_alleles_wless_than_number_of_supported_alleles(self):
        histograms = get_allele_histograms()
        noise_table = get_noise_table()
        alleles_1 = calculate_alleles(histograms[1], noise_table, 5)
        self.assertEqual(len(alleles_1.repeat_lengths), 1)

    def test_hard_support_threshold(self):
        histograms = get_allele_histograms()
        noise_table = get_noise_table()
        alleles_2 = calculate_alleles(histograms[2], noise_table, 5)
        self.assertEqual(len(alleles_2.repeat_lengths), 1)
        alleles_3 = calculate_alleles(histograms[3], noise_table, 5)
        self.assertEqual(len(alleles_3.repeat_lengths), 0)



        alleles_6 = calculate_alleles(histograms[5], noise_table, 5)
        self.assertEqual(int(alleles_6.log_likelihood), -67)
        alleles_5 = calculate_alleles(histograms[4], noise_table, 5)
        print(alleles_5)
        print(alleles_5.frequencies)
        print(histograms[4])
        self.assertEqual(len(alleles_5.repeat_lengths), 1)
        # alleles_7 = calculate_alleles(histograms[6], noise_table, 5)
        # sorted_freqs = sorted(alleles_7.frequencies)
        # print(sorted_freqs)
        # self.assertTrue(abs(sorted_freqs[1] - 0.230181 ) < .01)
        # self.assertTrue(abs(sorted_freqs[2] - 0.731404) < .01)
        # self.assertEqual(int(alleles_7.log_likelihood), -102)




if __name__ == '__main__':
    unittest.main()
