import os, unittest

from scipy.stats import friedmanchisquare

from src.Entry.SingleFileBatches import run_single_histogram
from tests.testing_utils.read_results import ResultsReader
from tests.testing_utils.sam_utils import FakeRead, create_MD_string, write_seq
from tests.testing_utils.self_contained_utils import run_msmutect_from_cmd, locus_file_path, sample_bams_path, \
    test_results_path, locus_file_path_strict, real_locus_file_path
from tests.testing_utils.write_accurate_sam_file import create_faithful_bam_16520, full_match_read, one_del_read, \
    one_insertion_read, one_insertion_2_deletion_read, wrong_sequence_insertion, deletion_over_end, \
    deletion_over_beginning, snp_pre_locus, snp_way_post_locus, deletion_with_SNP

j = os.path.join


class TestSamUtils(unittest.TestCase):

    def test_md_str_creation(self):
        reads = [
            FakeRead(1, "101M", "A"*101, []), # regular read
            FakeRead(1, "10M3D91M", "A" * 101, ),  # 3 base deletion
            FakeRead(1, "10M6D91M", "A"*101)
        ]

        modifications = [
            [],
            ["TTT"],
            ["TACTAC"],
        ]

        results = [
            "101",
            "10^TTT91",
            "10^TACTAC91"
        ]
        for i in range(len(reads)):
            read = reads[i]
            res = results[i]
            sub = modifications[i]
            md_str = create_MD_string(read, sub)
            self.assertEqual(md_str, "MD:Z:"+res)

    def test_seq_creation(self):
        reads_mods = [
            # (FakeRead(12, "10M3D91M"), ["TAC"]),
            (FakeRead(12, "5M1X3M3D92M", insertions_snps=["N"]), ["G", "CTA"]),
            (FakeRead(12, "5M6D96M"), ["ACT"*2]),
            (FakeRead(12, "20M3I2M6D79M", insertions_snps=["GGG"]), ["GGG", "TAC"*2])
        ]

        for r in reads_mods:
            read = r[0]
            seq, mods = write_seq(read.read_start, read.cigar_str, read.insertions_snps, base_sequence="ACT"*100, base_sequence_position=5)
            self.assertEqual(mods, r[1])


