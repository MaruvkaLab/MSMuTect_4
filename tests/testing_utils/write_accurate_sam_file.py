import re
import shutil, os
from dataclasses import dataclass
from typing import List

from tests.testing_utils.sam_utils import FakeRead, create_new_bam
from tests.testing_utils.sample_sequences import real_seq
from tests.testing_utils.self_contained_utils import sample_bams_path, header_only_sam


def full_match_read():
    return FakeRead(16_600, "101M")

def one_del_read():
    return FakeRead(16_600, "22M3D89M")

def one_insertion_read():
    return FakeRead(16_600, "20M3I86M", insertions_snps=["GCT"])

def one_insertion_2_deletion_read():
    return FakeRead(16_600, "20M3I2M6D86M", insertions_snps=["GCT"])

def wrong_sequence_insertion():
    return FakeRead(16_600, "20M3I88M", insertions_snps=["CGT"])

def deletion_over_end():
    return FakeRead(16_600, "31M3D70M")

def deletion_over_beginning():
    return FakeRead(16_600, "19M3D82M")

def snp_pre_locus():
    return FakeRead(16_600, "19M1X81M", insertions_snps=["N"])

def snp_post_locus():
    return FakeRead(16_600, "33M1X67M", insertions_snps=["N"])

def snp_way_post_locus():
    return FakeRead(16_600, "36M1X66M", insertions_snps=["N"])

def deletion_with_SNP():
    return FakeRead(16_600, "20M3D1X80M", insertions_snps=["N"])


def create_faithful_bam_16520(name: str, fake_reads: List[FakeRead]):
    create_new_bam(name, fake_reads, base_seq=real_seq(), base_seq_pos=16_520)

