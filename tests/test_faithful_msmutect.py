import os, unittest
from typing import List


from src.Entry.SingleFileBatches import run_single_histogram
from tests.testing_utils.read_results import ResultsReader
from tests.testing_utils.sam_utils import FakeRead
from tests.testing_utils.self_contained_utils import run_msmutect_from_cmd, locus_file_path, sample_bams_path, \
    test_results_path, locus_file_path_strict, real_locus_file_path, extended_locus_file_path
from tests.testing_utils.write_accurate_sam_file import create_faithful_bam_16520, full_match_read, one_del_read, \
    one_insertion_read, one_insertion_2_deletion_read, wrong_sequence_insertion, deletion_over_end, \
    deletion_over_beginning, snp_pre_locus, snp_way_post_locus, deletion_with_SNP, snp_post_locus, deletion_past_end
from collections import namedtuple

j = os.path.join

Case = namedtuple("Case", ["name", "reads", "lengths", "motif_support", "noisy"])


class TestFaithfulMSMuTect(unittest.TestCase):

    def run_case(self, name: str, fake_reads: List[FakeRead], locus_file_path=real_locus_file_path()):
        create_faithful_bam_16520(name, fake_reads)
        output = os.path.join(test_results_path(), name)
        cmd = f"-l {locus_file_path} -H -S {j(sample_bams_path(), f'{name}.bam')} -O {output} -f"
        print(cmd)
        run_msmutect_from_cmd(cmd)
        results_reader = ResultsReader(output + ".hist.tsv")
        first_line = next(results_reader)
        return first_line

    def test_cases(self):
        # first two fields are name and reads. Last 3 are the expected results
        cases = [
            Case('simple_faithful_test', [full_match_read()], [4], [1], False),
            Case('one_del', [one_del_read()], [3], [1], False),
            Case("one_insertion", [one_insertion_read()], [5], [1], False),
            Case("one_insertion_2_deletion", [one_insertion_2_deletion_read()], [3], [1], False),
            Case("wrong_insertion", [wrong_sequence_insertion()], [], [], False),
            Case("deletions_over_boundaries", [deletion_over_beginning(), deletion_over_end()], [], [], False),
            Case("noisy_snp_locus1", [snp_pre_locus() for _ in range(5)], [], [], True),
            Case("noisy_snp_locus2", [snp_post_locus() for _ in range(5)], [], [], True),
            Case("non_noisy_snp_locus", [snp_way_post_locus() for _ in range(5)], [4], [5], False),
            Case("deletion_with_snp", [deletion_with_SNP() for _ in range(4)], [], [], False)

        ]
        for case in cases:
            results = self.run_case(case.name, case.reads)
            self.assertEqual(case.lengths, results.motif_repeats, msg=f"{case.name}: Motif Lengths")
            self.assertEqual(case.motif_support, results.motif_repeat_support, msg=f"{case.name}: Motif Supports")
            self.assertEqual(case.noisy, results.noisy, msg=f"{case.name}: Noisy")

    def test_cases_fake_loci_file(self):
        # these cases use an extended locus so it can be impure
        # first two fields are name and reads. Last 3 are the expected results
        cases = [
            Case('deletion_past_end', [deletion_past_end()], [], [], False),
        ]
        for case in cases:
            results = self.run_case(case.name, case.reads, locus_file_path=extended_locus_file_path())
            self.assertEqual(case.lengths, results.motif_repeats, msg=f"{case.name}: Motif Lengths")
            self.assertEqual(case.motif_support, results.motif_repeat_support, msg=f"{case.name}: Motif Supports")
            self.assertEqual(case.noisy, results.noisy, msg=f"{case.name}: Noisy")










