from dataclasses import dataclass
from tests.E2E_tests.InputFile import InputFile
from tests.E2E_tests.load_config import local_normal_directory, local_tumor_directory, \
    normal_filename, docker_tumor_directory, tumor_filename, local_locus_directory, docker_locus_directory, \
    locus_filename, docker_output_directory, docker_normal_directory, local_output_directory, output_filename, \
    normal_hist_directory, normal_hist_filename, tumor_hist_directory, tumor_hist_filename


def default_normal_file() -> InputFile:
    return InputFile(local_file_directory=local_normal_directory(), docker_file_directory=docker_normal_directory(),
                     filename=normal_filename(), flag="-N")

def default_tumor_file() -> InputFile:
    return InputFile(local_file_directory=local_tumor_directory(), docker_file_directory=docker_tumor_directory(),
                     filename=tumor_filename(), flag="-T")

def default_locus_file() -> InputFile:
    return InputFile(local_file_directory=local_locus_directory(), docker_file_directory=docker_locus_directory(),
                     filename=locus_filename(), flag="-l")

def default_normal_hist_file() -> InputFile:
    return InputFile(local_file_directory=normal_hist_directory(), docker_file_directory=docker_normal_directory(),
                     filename=normal_hist_filename(), flag="-N")

def default_tumor_hist_file() -> InputFile:
    return InputFile(local_file_directory=tumor_hist_directory(), docker_file_directory=docker_tumor_directory(),
                     filename=tumor_hist_filename(), flag="-T")

def default_output_file() -> InputFile:
    return InputFile(local_file_directory=local_output_directory(), docker_file_directory=docker_output_directory(),
                     filename=output_filename(), flag="-O")

def default_single_file():
    ret = default_normal_file()
    ret.flag = "-S"
    return ret


@dataclass
class FileArguments:
    def __init__(self):
        self.normal_file: InputFile = default_normal_file()
        self.tumor_file: InputFile = default_tumor_file()
        self.locus_file: InputFile = default_locus_file()
        self.output_file: InputFile = default_output_file()
        self.single_file: InputFile = default_single_file()
        self.precompiled_normal_hist_file: InputFile = default_normal_hist_file()
        self.precompiled_tumor_hist_file: InputFile = default_tumor_hist_file()
