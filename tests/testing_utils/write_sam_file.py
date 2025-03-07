

from tests.testing_utils.sam_utils import FakeRead, create_new_bam
from tests.testing_utils.sample_sequences import seq_tri_repeat_full_purity
from tests.testing_utils.self_contained_utils import sample_bams_path, header_only_sam

# FakeRead = namedtuple("FakeRead", field_names=["read_start", "cigar_str"])


def main():

    snp_read = [FakeRead(9985, "15M1X6D85M", subsitutions=["N"]) for i in range(5)]
    create_new_bam("strict_test_1", [
        FakeRead(9985, "15M3D86M"),
        FakeRead(9985, "16M3D85M"),
        FakeRead(9985, "16M3I83M", subsitutions=["CTA"]),
        FakeRead(9985, "16M3I83M", subsitutions=["CAT"]),
        FakeRead(9985, "16M6D85M"),

    ]+snp_read, base_seq=seq_tri_repeat_full_purity())


    # create_new_bam("strict_test_2", [
    #     # these reads wont count for various reasons listed in the reads table
    #     FakeRead(9985, "14M3D87M"),
    #     FakeRead(9985, "70M15D31M"),
    #     FakeRead(9985, "13M90D88M"),
    #     FakeRead(9985, "12M1X88M", subsitutions=["N"]),
    #     FakeRead(9985, "40M10D71M"),
    #     FakeRead(9985, "40M5I71M", subsitutions=["AAAGA"]),
    #
    #     # these reads will count. see reads table
    #     FakeRead(9985, "90M1X10M", subsitutions=["N"]),
    #     FakeRead(9985, "72M2D29M", subsitutions=["N"])
    # ], create_seq_func=create_seq_mono_repeat_wone_impurity)




    # create_new_bam("test_locating_indels", [
    #     FakeRead(9985, "15M3I20M2D5M1D1I62M", subsitutions=["ACT", "T"]),
    #     FakeRead(9985, "14M2D87M"),
    #     FakeRead(9985, "74M2D27M"),
    #     FakeRead(9985, "15M60D86M"),
    #     FakeRead(9985, "10M70D91M"),
    #     FakeRead(9985, "1X1X87M1X1X10M", subsitutions=["N", "N", "N", "N"]),
    #
    #     FakeRead(9985, "101M") # just to make it easier to compare
    # ], create_seq_func=create_seq_high_entrophy)

    # fake_read_1 = [FakeRead(9985, "45M1X55M", subsitutions=["C"]) for i in range(3)] # snp at 10_031
    # fake_read_2 = [FakeRead(9986, "44M1X56M", subsitutions=["C"]) for i in range(2)] # should be single deletion in locus 3.]
    #
    # create_new_bam("test_edited_ref", fake_read_1+fake_read_2, create_seq_func=create_seq_high_entrophy)

    # fake_read_1 = [FakeRead(9985, "41M1X1M1X57M", subsitutions=["A", "A"]) for i in range(5)] # should be single deletion in locus 3.]
    # create_new_bam("test_elaborate_ref_based",
    #     fake_read_1+ \
    #     [
    #     FakeRead(9985, "66M3D35M"), # RD 2, should be single deletion in locus 4
    #     FakeRead(9985, "62M3D20M4D19M"), # RD 3, should be a deletion in locus 3 (and not 4)
    #     FakeRead(9985, "41M1X1M1X5M4D53M", subsitutions=["A", "A"]), # RD 4 should be a deletion of 1 of locus 1. 48M4D53M
    #     FakeRead(9985, "50M20D51M"), # RD 5 should wipe out entire third locus
    #     FakeRead(9985, "14M10D87M"), # RD 6, # should delete 8 bases from last locus
    #     FakeRead(9985, "70M20D31M"), # RD 7, should delete 6 bases from last locus
    #     FakeRead(9992, "42M3D59M"), # RD 8, double deletion in locus 2
    #     FakeRead(9992, "101M"), # RD 9, reference, does not support locus 1 (flanking)
    # ])

    # create_new_bam("dentist",
    #     [
    #     FakeRead(10_000, "25M3I25M1X5D47M", subsitutions=["TTT", "A"]),
    #     FakeRead(10_000, "50M20D51M"),
    #     FakeRead(10_000, "21M50D79M"),
    #     FakeRead(10_000, "40M2I59M", subsitutions=["TT"]),
    #     FakeRead(10_000, "101M"),
    #     ], create_seq_func=create_seq_high_entrophy)

    # create_new_bam("test_char_count_annotated", [
    #         FakeRead(1, "101M", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
    #         FakeRead(10, "13M4D88M"),
    #         FakeRead(15, "12M3D89M"),
    #         FakeRead(15, "10M3D91M"),
    #         FakeRead(15, "10M3I88M"),
    #         FakeRead(15, "20M10D81M"),
    #         FakeRead(15, "19M10D82M"),
    #         FakeRead(15, "10M15D91M"),
    #
    # ])
    # create_new_bam("test_annotated_locus", [
    #     FakeRead(1, "101M",
    #              "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), # shouldnt be used for consensus: total match
    #     FakeRead(1, "26M2D2S73M",
    #              "AAAAAAAAAAAAAAAAAAAAAAAAAAGGACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA"),# shouldnt be used for consensus: deletion
    #     FakeRead(1, "26M2I73M",
    #              "AAAAAAAAAAAAAAAAAAAAAAAAAAGGACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA"), # shouldnt be used for consensus: insertion
    #     FakeRead(1, "26M2X73M",
    #              "AAAAAAAAAAAAAAAAAAAAAAAAAAGGACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA"),
    #     FakeRead(1, "26M2X73M",
    #              "AAAAAAAAAAAAAAAAAAAAAAAAAAGGACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA"),
    #     FakeRead(1, "26M2X73M",
    #              "AAAAAAAAAAAAAAAAAAAAAAAAAAGGACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA"),
    #     FakeRead(1, "26M2X73M",
    #              "AAAAAAAAAAAAAAAAAAAAAAAAAAGGACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA"),
    #     FakeRead(1, "26M2X73M",
    #              "AAAAAAAAAAAAAAAAAAAAAAAAAAGGACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA"),
    #     FakeRead(1, "26M2X73M",
    #              "AAAAAAAAAAAAAAAAAAAAAAAAAAGGACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA"),
    #     FakeRead(1, "26M2X73M",
    #              "AAAAAAAAAAAAAAAAAAAAAAAAAAGGACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA"),
    #
    # ])

    # create_new_bam("test_extract_locus", [
    #     FakeRead(1, "101M",
    #              "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA"),
    #     FakeRead(6, "101M",
    #              "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA"),
    #     FakeRead(6, "12M2D89M",
    #              "ACGTACGTACGTGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG"),
    # ])

    # create_new_bam("test_full_pipe", [
    #     FakeRead(9985, "101M",
    #              "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA"),
    #     FakeRead(9987, "101M",
    #              "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA"),
    #     FakeRead(9988, "101M",
    #              "ACGTACGTACGTGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG"),
    #
    # ])

    # create_new_bam("test_full_pipe_complex",[
    #     FakeRead(9989, "101M",
    #              "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACACACACACAAAAACGACGACGACGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
    #     FakeRead(9990, "101M",
    #              "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACACACACACAAAAACGACGACGACGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
    #     FakeRead(9990, "101M",
    #              "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACACACACACAAAAACGACGACGACGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
    #     FakeRead(9990, "101M",
    #              "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACACACACACAAAAACGACGACGACGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
    #     FakeRead(9990, "101M",
    #              "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACACACACACAAAAACGACGACGACGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
    #     FakeRead(9990, "36M2D65M",
    #              "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACACACACAAAAACGACGACGACGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), # deletion should get assigned to first locus
    # ]) # so it does a count for the


    # create_new_bam("mapping_small_locus", [
    #         FakeRead(9_964, "101M"),  # shouldn't map (misses on flanking)
    #         FakeRead(9_965, "101M"),  # should map
    #         FakeRead(10_035, "101M"),  # should map
    #         FakeRead(10_036, "101M"),  # shouldn't map
    # ])
    #
    #
    # create_new_bam("indels", [
    #         FakeRead(10_000, "35M30D66M"),  # knock out entire smallest locus. repeat_length=0
    #
    #         FakeRead(10_035, "18M2D83M"),  # knock out end of smallest locus, repeat_length=9
    #         FakeRead(10_035, "19M10D82M"),  # knock out end and then some of smallest locus, repeat_length=9
    #         FakeRead(10_035, "19M3D82M"),  # knock out end and then some of smallest locus, repeat_length=9
    #         FakeRead(10_035, "9M3D92M"),  # knock out start of smallest locus, repeat_length=9
    #
    #         FakeRead(10_035, "10M1D91M"),  # knock out start of smallest locus, repeat_length=10
    #         FakeRead(10_035, "10M1D91M"),  # knock out start of smallest locus, repeat_length=10
    #
    #
    #         FakeRead(10_035, "10M1I90M"),  # insertion at start of smallest locus, repeat_length=12
    #         FakeRead(10_035, "20M1I80M"),  # insertion at end of smallest locus, repeat_length=12
    #         FakeRead(10_035, "15M1I85M"),  # insertion of middle of smallest locus, repeat_length=12
    #
    #
    #         FakeRead(10_035, "5M1I95M"),  # insertion before start of smallest locus, repeat_length=11
    #         FakeRead(10_035, "40M2D61M"),  # deletion after end of smallest locus, repeat_length=11
    #         FakeRead(10_035, "101M"),   # standard: rl=11
    #         FakeRead(10_035, "101M"),  # standard: rl=11
    #         FakeRead(10_035, "101M"),   # standard: rl=11
    #         FakeRead(10_035, "101M"),   # standard: rl=11
    #
    #     # 11_6, 9_4, 12_3, 10_2, 0_1
    # ])
    #
    #
    # create_new_bam_custom_length("multimapping_loci", [
    #         FakeRead(9990, "150M"),  # should align to all
    #         FakeRead(9991, "150M"), # should align to all except biggest
    #         FakeRead(10_001, "150M"), # should align to all except biggest
    #         FakeRead(10_026, "150M"), # should align to first and third
    #         FakeRead(10_040, "150M") # should not align to any
    # ], 150)


if __name__ == '__main__':
    main()