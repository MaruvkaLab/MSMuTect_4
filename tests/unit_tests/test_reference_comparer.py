import unittest, time

from src.GenomicUtils.Mutation import Mutation
from src.GenomicUtils.reference_locus_comparer import microsatellite_indel, single_ms_indel_determination



class TestAnnotatedLocus(unittest.TestCase):

    def test_single_ms_pattern_indel(self):
        self.assertTrue(single_ms_indel_determination("A", "A"))
        self.assertTrue(single_ms_indel_determination("ACT", "ACT"))
        self.assertTrue(single_ms_indel_determination("CAC", "CCA"))
        self.assertTrue(single_ms_indel_determination("CAC", "ACC"))

        self.assertFalse(single_ms_indel_determination("ACT", "ATC"))
        self.assertFalse(single_ms_indel_determination("CCT", "ATC"))


    def test_ms_indel(self):
        mono_repeat_1 = Mutation("A", insertion=True)
        self.assertEqual(microsatellite_indel(mono_repeat_1, "A"), 1)
        #
        mono_repeat_3 = Mutation("AAA", deletion=True)
        self.assertEqual(microsatellite_indel(mono_repeat_3, "A"), -3)
        #
        quad_repeat_2 = Mutation("ACTGACTG", insertion=True)
        self.assertEqual(microsatellite_indel(quad_repeat_2, "ACTG"), 2)
        #
        tri_repeat_complex = Mutation("CCACCACCACCA", insertion=True)
        self.assertEqual(microsatellite_indel(tri_repeat_complex, "CAC"), 4)
        #
        problematic_incomplete = Mutation("ACTGACT", insertion=True)
        self.assertEqual(microsatellite_indel(problematic_incomplete, "ACTG"), 0)

        problematic_incomplete_2 = Mutation("AACAA", deletion=True)
        self.assertEqual(microsatellite_indel(problematic_incomplete_2, "A"), 0)

        problematic_incomplete_2 = Mutation("TCA", deletion=True)
        self.assertEqual(microsatellite_indel(problematic_incomplete_2, "ATC"), -1)




if __name__ == '__main__':
    unittest.main()
