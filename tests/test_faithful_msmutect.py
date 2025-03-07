import os, unittest
from typing import List

from scipy.stats import friedmanchisquare

from src.Entry.SingleFileBatches import run_single_histogram
from tests.testing_utils.read_results import ResultsReader
from tests.testing_utils.sam_utils import FakeRead
from tests.testing_utils.self_contained_utils import run_msmutect_from_cmd, locus_file_path, sample_bams_path, \
    test_results_path, locus_file_path_strict, real_locus_file_path
from tests.testing_utils.write_accurate_sam_file import create_faithful_bam_16520, full_match_read, one_del_read, \
    one_insertion_read, one_insertion_2_deletion_read, wrong_sequence_insertion, deletion_over_end, \
    deletion_over_beginning, snp_pre_locus, snp_way_post_locus, deletion_with_SNP
from collections import namedtuple

j = os.path.join


class TestFaithfulMSMuTect(unittest.TestCase):

    def run_case(self, name: str, fake_reads: List[FakeRead]):
        create_faithful_bam_16520(name, fake_reads)
        output = os.path.join(test_results_path(), name)
        cmd = f"-l {real_locus_file_path()} -H -S {j(sample_bams_path(), f'{name}.bam')} -O {output} -f"
        print(cmd)
        run_msmutect_from_cmd(cmd)


    def test_cases(self):
        # first two fields are name and reads. Last 3 are the expected results
        Case = namedtuple("Case", ["name", "reads", "lengths", "motif_support", "noisy"])
        cases = [
            Case('simple_faithful_test', [full_match_read()], [4], [1], False)
        ]
        for case in cases:
            self.run_case(case.name, case.reads)
            self.assertEqual()


    def test_can_run(self):
        # one read supporting the regular locus
        name = 'simple_faithful_test'
        create_faithful_bam_16520(name, [full_match_read()])
        output = os.path.join(test_results_path(), name)
        cmd = f"-l {real_locus_file_path()} -H -S {j(sample_bams_path(), f'{name}.bam')} -O {output} -f"
        print(cmd)
        run_msmutect_from_cmd(cmd)

    def test_one_del(self):
        # one read supporting single motif deletion
        name = 'one_del'
        create_faithful_bam_16520(name, [one_del_read()])
        output = os.path.join(test_results_path(), name)
        cmd = f"-l {real_locus_file_path()} -H -S {j(sample_bams_path(), f'{name}.bam')} -O {output} -f"
        print(cmd)
        run_msmutect_from_cmd(cmd)
        results_reader = ResultsReader(output+".hist.tsv")
        first_line = next(results_reader)
        self.assertEqual(len(first_line.motif_repeats), 1)
        self.assertEqual(first_line.motif_repeats[0], 3)


    def test_one_insertion(self):
        # one read supporting single motif deletion
        name = "one_insertion"
        create_faithful_bam_16520(name, [one_insertion_read()])
        output = os.path.join(test_results_path(), name)
        cmd = f"-l {real_locus_file_path()} -H -S {j(sample_bams_path(), f'{name}.bam')} -O {output} -f"
        print(cmd)
        run_msmutect_from_cmd(cmd)
        results_reader = ResultsReader(output+".hist.tsv")
        first_line = next(results_reader)
        self.assertEqual(len(first_line.motif_repeats), 1)
        self.assertEqual(first_line.motif_repeats[0], 5)

    def test_1_insertion_2_deletion(self):
        # one read supporting single motif deletion
        name = "one_insertion_2_deletion"
        create_faithful_bam_16520(name, [one_insertion_2_deletion_read()])
        output = os.path.join(test_results_path(), name)
        cmd = f"-l {real_locus_file_path()} -H -S {j(sample_bams_path(), f'{name}.bam')} -O {output} -f"
        print(cmd)
        run_msmutect_from_cmd(cmd)
        results_reader = ResultsReader(output+".hist.tsv")
        first_line = next(results_reader)
        self.assertEqual(len(first_line.motif_repeats), 1)
        self.assertEqual(first_line.motif_repeats[0], 3)

    def test_wrong_insertion(self):
        # one read supporting single motif deletion
        name = "wrong_insertion"
        create_faithful_bam_16520(name, [wrong_sequence_insertion()])
        output = os.path.join(test_results_path(), name)
        cmd = f"-l {real_locus_file_path()} -H -S {j(sample_bams_path(), f'{name}.bam')} -O {output} -f"
        print(cmd)
        run_msmutect_from_cmd(cmd)
        results_reader = ResultsReader(output+".hist.tsv")
        first_line = next(results_reader)
        self.assertEqual(len(first_line.motif_repeats), 0)

    def test_deletions_over_beginning_and_end(self):
        # one read supporting single motif deletion
        name = "deletions_over_boundaries"
        create_faithful_bam_16520(name, [deletion_over_beginning(), deletion_over_end()])
        output = os.path.join(test_results_path(), name)
        cmd = f"-l {real_locus_file_path()} -H -S {j(sample_bams_path(), f'{name}.bam')} -O {output} -f"
        print(cmd)
        run_msmutect_from_cmd(cmd)
        results_reader = ResultsReader(output+".hist.tsv")
        first_line = next(results_reader)
        self.assertEqual(len(first_line.motif_repeats), 0)

    def test_noisy_snp_locus2(self):
        # one read supporting single motif deletion
        name = "noisy_snp_locus2"
        create_faithful_bam_16520(name, [snp_pre_locus() for _ in range(5)])
        output = os.path.join(test_results_path(), name)
        cmd = f"-l {real_locus_file_path()} -H -S {j(sample_bams_path(), f'{name}.bam')} -O {output} -f"
        print(cmd)
        run_msmutect_from_cmd(cmd)
        results_reader = ResultsReader(output+".hist.tsv")
        first_line = next(results_reader)
        self.assertEqual(len(first_line.motif_repeats), 0)
        self.assertTrue(first_line.noisy)

    def test_noisy_snp_locus(self):
        # one read supporting single motif deletion
        name = "noisy_snp_locus"
        create_faithful_bam_16520(name, [snp_pre_locus() for _ in range(5)])
        output = os.path.join(test_results_path(), name)
        cmd = f"-l {real_locus_file_path()} -H -S {j(sample_bams_path(), f'{name}.bam')} -O {output} -f"
        print(cmd)
        run_msmutect_from_cmd(cmd)
        results_reader = ResultsReader(output+".hist.tsv")
        first_line = next(results_reader)
        self.assertEqual(len(first_line.motif_repeats), 0)
        self.assertTrue(first_line.noisy)

    def test_non_noisy(self):
        # one read supporting single motif deletion
        name = "non_noisy_snp_locus"
        create_faithful_bam_16520(name, [snp_way_post_locus() for _ in range(5)])
        output = os.path.join(test_results_path(), name)
        cmd = f"-l {real_locus_file_path()} -H -S {j(sample_bams_path(), f'{name}.bam')} -O {output} -f"
        print(cmd)
        run_msmutect_from_cmd(cmd)
        results_reader = ResultsReader(output+".hist.tsv")
        first_line = next(results_reader)
        self.assertEqual(len(first_line.motif_repeats), 1)
        self.assertEqual(first_line.motif_repeats[0], 4)
        self.assertFalse(first_line.noisy)

    # def test_deletion_with_snp(self):
    #     # one read supporting single motif deletion
    #     name = "deletion_with_snp"
    #     create_faithful_bam_16520(name, [deletion_with_SNP()])
    #     output = os.path.join(test_results_path(), name)
    #     cmd = f"-l {real_locus_file_path()} -H -S {j(sample_bams_path(), f'{name}.bam')} -O {output} -f"
    #     print(cmd)
    #     run_msmutect_from_cmd(cmd)
    #     results_reader = ResultsReader(output+".hist.tsv")
    #     first_line = next(results_reader)
    #     self.assertEqual(len(first_line.motif_repeats), 0)
    #     self.assertFalse(first_line.noisy)








