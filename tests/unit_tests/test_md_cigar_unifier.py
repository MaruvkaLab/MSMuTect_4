import os, unittest

from src.GenomicUtils.Mutation import Mutation
from src.GenomicUtils.md_cigar_parser import split_md_string_into_tokens, split_MD_into_ops, MD_OP, MD_TUPLE


class TestStrictMSMuTect(unittest.TestCase):

    def test_split_md_string(self):
        # self.assertListEqual([(1,2)], [(1,2)])
        self.assertListEqual(split_md_string_into_tokens("103A^CCA13A"), ["103", "A", "^", "C", "C", "A", "13", "A"])
        self.assertListEqual(split_md_string_into_tokens("40^TTTA3A107"), ["40", "^", "T", "T", "T", "A", "3", "A", "107"])


        self.assertListEqual(split_MD_into_ops("103A^CCA13A"), [MD_TUPLE(MD_OP.MATCH, 103, "N"),
                                                                    MD_TUPLE(MD_OP.SUBSITUTION, 1, "N"),
                                                                MD_TUPLE(MD_OP.DELETION, 3, "CCA"),
                                                                MD_TUPLE(MD_OP.MATCH, 13, "N"),
                                                                MD_TUPLE(MD_OP.SUBSITUTION, 1, "N"),
                                                                ])

        self.assertListEqual(split_MD_into_ops("110^CC"), [MD_TUPLE(MD_OP.MATCH, 110, "N"),
                                                           MD_TUPLE(MD_OP.DELETION, 2, "CC")])

        self.assertListEqual(split_MD_into_ops("40^TTTA3A107"), [
                                                                MD_TUPLE(MD_OP.MATCH, 40, "N"),
                                                                MD_TUPLE(MD_OP.DELETION, 4, "TTTA"),
                                                                MD_TUPLE(MD_OP.MATCH, 3, "N"),
                                                                MD_TUPLE(MD_OP.SUBSITUTION, 1, "N"),
                                                                MD_TUPLE(MD_OP.MATCH, 107, "N")
                                                                ])

