import os, json
from typing import Dict


def current_directory()->str:
    return os.path.dirname(os.path.realpath(__file__))

def config_path()->str:
    return os.path.join(current_directory(), "config.json")

def load_config()->Dict[str, str]:
    with open(config_path()) as f:
        config = json.load(f)
    return config

def output_directory() -> str:
    return os.path.join(current_directory(), "output")

def get_parameter(parameter_name: str) -> str:
    config = load_config()
    ret = config[parameter_name]
    return ret

def directory_of_parameter(parameter_name: str) -> str:
    parameter = get_parameter(parameter_name)
    return os.path.dirname(parameter)

def filename_of_parameter(parameter_name: str) -> str:
    parameter = get_parameter(parameter_name)
    return os.path.basename(parameter)

def local_locus_directory()->str:
    return directory_of_parameter("locus_file")

def local_tumor_directory()->str:
    return directory_of_parameter("full_tumor_bam")

def local_normal_directory()->str:
    return directory_of_parameter("full_normal_bam")

def local_output_directory() -> str:
    return directory_of_parameter("local_output")

def docker_locus_directory()->str:
    return get_parameter("docker_locus_directory")

def docker_tumor_directory()->str:
    return get_parameter("docker_tumor_directory")

def docker_normal_directory()->str:
    return get_parameter("docker_normal_directory")

def docker_output_directory()->str:
    return get_parameter("docker_output_directory")

def tumor_filename()->str:
    return filename_of_parameter("full_tumor_bam")

def normal_filename()->str:
    return filename_of_parameter("full_normal_bam")

def locus_filename()->str:
    return filename_of_parameter("locus_file")

def output_filename() -> str:
    return filename_of_parameter("local_output")

def normal_intermediate_file()->str:
    raise NotImplementedError

def tumor_intermediate_file()->str:
    raise NotImplementedError

def executable_path()->str:
    return get_parameter("executable_path")

def bash_path()->str:
    return get_parameter("bash_script_path")

def normal_hist_directory()->str:
    return directory_of_parameter("normal_hist_tsv")

def tumor_hist_directory()->str:
    return directory_of_parameter("tumor_hist_tsv")

def normal_hist_filename():
    return filename_of_parameter("normal_hist_tsv")

def tumor_hist_filename()->str:
    return filename_of_parameter("tumor_hist_tsv")

def ground_truth_hist_file() -> str:
    return os.path.join(normal_hist_directory(), normal_hist_filename())

if __name__ == "__main__":
    print(get_parameter("executable_path"))