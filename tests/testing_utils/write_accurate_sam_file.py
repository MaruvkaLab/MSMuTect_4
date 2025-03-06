import re
import shutil, os
from dataclasses import dataclass
from typing import List

from tests.testing_utils.sam_utils import FakeRead
from tests.testing_utils.self_contained_utils import sample_bams_path, header_only_sam




def real_seq():
    #10,001-10,468
    return "TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAAACCCTAAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCTAACCCCTAACCCTAACCCTAACCCTACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCT"


def main():
    snp_read = [FakeRead(9985, "15M1X6D85M", subsitutions=["N"]) for i in range(5)]
    create_new_bam("strict_test_1", [
        FakeRead(9985, "15M3D86M"),
        FakeRead(9985, "16M3D85M"),
        FakeRead(9985, "16M3I83M", subsitutions=["CTA"]),
        FakeRead(9985, "16M3I83M", subsitutions=["CAT"]),
        FakeRead(9985, "16M6D85M"),

    ]+snp_read, create_seq_func=create_seq_tri_repeat_full_purity)