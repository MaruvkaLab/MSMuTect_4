import unittest, time
from tests.testing_utils.self_contained_utils import locus_file_path
from src.GenomicUtils.LocusFile import LociManager


class TestLociManager(unittest.TestCase):

    # def test_locus_load(self):
    #     init_start = time.process_time()
    #     manager = LociManager("/home/avraham/MaruvkaLab/msmutect_runs/data/hg19_all_loci_sorted_new", 1_000_000)
    #     self.assertLess(time.process_time()-init_start, 4)
    #     init_start = time.process_time()
    #     _ = manager.get_batch(1_000_000)
    #     self.assertLess(time.process_time()-init_start, 8)

    def test_load_single(self):
        manager = LociManager(locus_file_path(), 0)
        l = manager.get_batch(1)[0]
        self.assertTrue(True)

    def test_load_seung_won(self):
        manager = LociManager("/home/avraham/MaruvkaLab/msmutect_runs/data/seung_won.tsv", 0)
        l = manager.get_batch(1)[0]
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
