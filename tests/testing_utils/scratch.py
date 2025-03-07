from tests.testing_utils.write_accurate_sam_file import create_faithful_bam_16520, one_insertion_2_deletion_read

name = "croc"
create_faithful_bam_16520(name, [one_insertion_2_deletion_read()])