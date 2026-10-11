import os

from pysam import AlignmentFile

from src.GenomicUtils.md_cigar_parser import extract_md_str
from tests.testing_utils.self_contained_utils import sample_bams_path


def all_reads_from_bam_file(bam_file_fp: str):
    BAM_handle = AlignmentFile(bam_file_fp, "rb")
    a=BAM_handle.fetch("1", start=1, multiple_iterators=False)
    all_reads = []
    while True:
        try:
            all_reads.append(next(a))
        except StopIteration:
            return all_reads

def main():
    # a = all_reads_from_bam_file("/home/avraham/MaruvkaLab/msmutect_runs/HordeFeatureUpdate/1kreads.bam")
    a = all_reads_from_bam_file(os.path.join(sample_bams_path(), 'strict_test_2.bam'))
    for j in a:
        print(extract_md_str(j))
        # if j.is_reverse:
        #     print("reverse")
        #     print(j.seq)
        #     # print(j.get_aligned_pairs())
        #     print("***************************")
        # else:
        #     print("croc")
        #     print(j.seq)
        #     print("***************************")

if __name__ == '__main__':
    main()