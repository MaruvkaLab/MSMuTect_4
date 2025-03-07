import os.path


def sample_bams_path():
    p = os.path.dirname
    j = os.path.join
    return j(p(p(__file__)), "sample_bams")


def header_only_sam():
    return os.path.join(sample_bams_path(), "header_only.sam")


def locus_file_path():
    return os.path.join(sample_bams_path(), "fake_sample_locus_sorted.tsv")

def real_locus_file_path():
    return os.path.join(sample_bams_path(), "real_locus.phobos")

def extended_locus_file_path():
    return os.path.join(sample_bams_path(), "real_locus_impure.phobos")

def locus_file_path_strict():
    return os.path.join(sample_bams_path(), "fake_loci_strict.tsv")


def msmutect_path():
    p = os.path.dirname
    return p(p(p(__file__)))


def msmutect_executable_path():
    return os.path.join(msmutect_path(), "msmutect.sh")


def run_msmutect_from_cmd(args: str):
    ret_code = os.system(msmutect_executable_path() + f" {args}")
    if ret_code != 0:
        raise RuntimeError("MSMuTect crashed")


def test_results_path():
    return os.path.join(msmutect_path(), "tests", "test_results")