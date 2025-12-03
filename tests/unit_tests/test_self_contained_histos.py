import os, unittest

from tests.testing_utils.generate_histograms import histogram_histograms
from tests.testing_utils.read_results import ResultsReader
from tests.testing_utils.self_contained_utils import run_msmutect_from_cmd, locus_file_path, sample_bams_path, \
    test_results_path


class CrocTrap:#TestHistogram(unittest.TestCase):

    def test_mapping(self):
        j = os.path.join
        mapping_small_locus_regular_flanking = os.path.join(test_results_path(), 'mapping')
        run_msmutect_from_cmd(f"-l {locus_file_path()} -H -S {j(sample_bams_path(), 'mapping_small_locus.bam')}"
                              f" -O {mapping_small_locus_regular_flanking} -f")
        results_reader = ResultsReader(mapping_small_locus_regular_flanking+".hist.tsv")
        first_line = next(results_reader)
        self.assertEqual(len(first_line.motif_repeat_support), 1)
        self.assertEqual(first_line.motif_repeat_support[0], 2)

    def test_mapping_custom_flanking(self):
        j = os.path.join
        mapping_small_locus_regular_flanking = os.path.join(test_results_path(), 'mapping_custom_flanking')
        run_msmutect_from_cmd(f"-l {locus_file_path()} -H -S {j(sample_bams_path(), 'mapping_small_locus.bam')}"
                              f" -O {mapping_small_locus_regular_flanking} -f --flanking 9")
        results_reader = ResultsReader(mapping_small_locus_regular_flanking + ".hist.tsv")
        first_line = next(results_reader)
        self.assertEqual(len(first_line.motif_repeat_support), 1)
        self.assertEqual(first_line.motif_repeat_support[0], 4)

    def test_rounding_winteger_only_flag(self):
        rounded_histogram = histogram_histograms()[0]
        self.assertEqual(0, rounded_histogram.rounded_repeat_lengths[11])
        self.assertEqual(0, rounded_histogram.rounded_repeat_lengths[12])

        rounded_histogram_1 = histogram_histograms()[1]
        self.assertEqual(5, rounded_histogram_1.rounded_repeat_lengths[11])
        self.assertEqual(3, rounded_histogram_1.rounded_repeat_lengths[12])

        rounded_histogram_2 = histogram_histograms()[2]
        self.assertEqual(3, rounded_histogram_2.rounded_repeat_lengths[5])

    def test_indels(self):
        j = os.path.join
        indel_results = os.path.join(test_results_path(), 'indels')
        run_msmutect_from_cmd(f"-l {locus_file_path()} -H -S {j(sample_bams_path(), 'indels.bam')}"
                              f" -O {indel_results} -f")
        results_reader = ResultsReader(indel_results + ".hist.tsv")
        first_line = next(results_reader)
        self.assertEqual(len(first_line.motif_repeat_support), 5)

        # 11_6, 9_4, 12_3, 10_2, 0_1

        self.assertEqual(first_line.motif_repeats[0], 11)
        self.assertEqual(first_line.motif_repeats[1], 9)
        self.assertEqual(first_line.motif_repeats[2], 12)
        self.assertEqual(first_line.motif_repeats[3], 10)
        self.assertEqual(first_line.motif_repeats[4], 0)

        self.assertEqual(first_line.motif_repeat_support[0], 6)
        self.assertEqual(first_line.motif_repeat_support[1], 4)
        self.assertEqual(first_line.motif_repeat_support[2], 3)
        self.assertEqual(first_line.motif_repeat_support[3], 2)
        self.assertEqual(first_line.motif_repeat_support[4], 1)

        second_line = next(results_reader)
        self.assertEqual(second_line.motif_repeats[0], 0)
        self.assertEqual(second_line.motif_repeat_support[0], 1)

        third_line = next(results_reader)
        self.assertEqual(third_line.motif_repeats[0], 0)
        self.assertEqual(third_line.motif_repeat_support[0], 1)

    def test_multimapping_loci(self):
        j = os.path.join
        multimapping_results = os.path.join(test_results_path(), 'multimapping')
        run_msmutect_from_cmd(f"-l {locus_file_path()} -H -S {j(sample_bams_path(), 'multimapping_loci.bam')}"
                              f" -O {multimapping_results} -f")
        results_reader = ResultsReader(multimapping_results + ".hist.tsv")

        mapping_reads = [4, 3, 4, 1]
        for i in range(4):
            current_line = next(results_reader)
            self.assertEqual(mapping_reads[i], current_line.motif_repeat_support[0], f"Failed on {i}")


if __name__ == '__main__':
    unittest.main()

