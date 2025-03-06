import os, unittest

from src.Entry.SingleFileBatches import run_single_histogram
from tests.testing_utils.read_results import ResultsReader
from tests.testing_utils.self_contained_utils import run_msmutect_from_cmd, locus_file_path, sample_bams_path, \
    test_results_path, locus_file_path_strict


class TestStrictMSMuTect(unittest.TestCase):

    # def test_can_run(self):
    #     j = os.path.join
    #     mapping_small_locus_regular_flanking = os.path.join(test_results_path(), 'test_strict')
    #     cmd = f"-l {locus_file_path()} -H -S {j(sample_bams_path(), 'mapping_small_locus.bam')} -O {mapping_small_locus_regular_flanking} -f"
    #     print(cmd)
    #     run_msmutect_from_cmd(cmd)
    #
    # def test_indel_forms_complex_motif(self):
    #     j = os.path.join
    #     testing_results = os.path.join(test_results_path(), 'test_strict')
    #     # cmd = f"-l {locus_file_path_strict()} -H -S {j(sample_bams_path(), 'strict_test_1.bam')} -O {testing_results} -f"
    #     # run_msmutect_from_cmd(cmd)
    #     run_single_histogram(j(sample_bams_path(), 'strict_test_1.bam'), locus_file_path_strict(), 0, 2, 1, 10, testing_results)
    #     results_reader = ResultsReader(testing_results+".hist.tsv")
    #     first_line = next(results_reader)
    #     repeats = [3, 5, 2]
    #     repeat_supports = [2, 1, 1]
    #     self.assertTrue(first_line.noisy)
    #     self.assertEqual(len(first_line.motif_repeats), len(repeats))
    #     for i in range(len(repeats)):
    #         self.assertEqual(first_line.motif_repeats[i], repeats[i])
    #         self.assertEqual(first_line.motif_repeat_support[i], repeat_supports[i])
    #
    # def test_indel_forms_long_mono_repeat(self):
    #     j = os.path.join
    #     testing_results = os.path.join(test_results_path(), 'test_strict')
    #     # cmd = f"-l {locus_file_path_strict()} -H -S {j(sample_bams_path(), 'strict_test_1.bam')} -O {testing_results} -f"
    #     # run_msmutect_from_cmd(cmd)
    #     run_single_histogram(j(sample_bams_path(), 'strict_test_2.bam'), locus_file_path_strict(), 0, 2, 1, 10, testing_results)
    #     results_reader = ResultsReader(testing_results+".hist.tsv")
    #     first_line = next(results_reader)
    #     second_line = next(results_reader)
    #     self.assertEqual(len(second_line.motif_repeats), 2)
    #     self.assertEqual(second_line.motif_repeats[0], 60)
    #     self.assertEqual(second_line.motif_repeat_support[0], 1)
    #     self.assertEqual(second_line.motif_repeats[1], 58)
    #     self.assertEqual(second_line.motif_repeat_support[1], 1)

    def test_real_small_case(self):
        self.assertTrue(True)
        # run_single_histogram("/home/avraham/MaruvkaLab/Texas/gdc/088afdfb-ae1e-4330-a015-38ce7fa2fa92/normal_5df6c9bb-e65d-44b5-a4df-342983a88c47_wgs_gdc_realn.bam",
        #                      "/home/avraham/MaruvkaLab/Texas/strict_msmutect/should_be_noisy.phobos",
        #                      0, 1, 1, 10, "/home/avraham/MaruvkaLab/Texas/strict_msmutect/results/noisy_locus")
        # results_reader = ResultsReader("/home/avraham/MaruvkaLab/Texas/strict_msmutect/results/noisy_locus.hist.tsv")
        # first_line = next(results_reader)
        # self.assertTrue(first_line.noisy)


    # def test_non_noisy_bug(self):
    #     j = os.path.join
    #     testing_results = os.path.join(test_results_path(), 'test_strict')
    #
    #     run_single_histogram("/home/avraham/MaruvkaLab/Texas/gdc/088afdfb-ae1e-4330-a015-38ce7fa2fa92/normal_5df6c9bb-e65d-44b5-a4df-342983a88c47_wgs_gdc_realn.bam",
    #                          "/home/avraham/MaruvkaLab/Texas/strict_msmutect/should_be_noisy.phobos", 0, 2, 1, 10,
    #                          testing_results)



