import os
from typing import List
from tests.E2E_tests.Enums import STAGE, EXECUTION_TYPE
from tests.E2E_tests.FileArguments import FileArguments
from tests.E2E_tests.InputFile import InputFile
from tests.E2E_tests.load_config import executable_path, bash_path, local_output_directory, output_filename



class CMDArguments:
    def __init__(self, executor: EXECUTION_TYPE, start_point: STAGE, end_point: STAGE, num_loci: int = None, num_cores: int = 1):
        file_arguments = FileArguments()
        input_files, run_flags  = self.extract_run_command(start_point, end_point, file_arguments)
        auxiliary_arguments = self.extract_auxiliary_arguments(num_loci, num_cores)

        start_commands = {
            EXECUTION_TYPE.EXECUTABLE: executable_path() + " ",
            EXECUTION_TYPE.BASH: bash_path() + " ",
            EXECUTION_TYPE.DOCKER: f"docker run {self.make_binding_arguments(input_files)} --rm msmutect-docker "
        }

        suffixes = {
            STAGE.HISTOGRAMS: ".hist.tsv",
            STAGE.ALLELES: ".all.tsv",
            STAGE.MUTATIONS_FULL: ".full.mut.tsv",
            STAGE.MUTATIONS_EFFICIENT: ".partial.mut.tsv",
        }

        output_suffix = suffixes[end_point]
        exe_cmd = start_commands[executor]
        self.command = exe_cmd + auxiliary_arguments + run_flags + self.format_input_files(executor, input_files) + " -f"
        self.local_output_file = os.path.join(local_output_directory(), file_arguments.output_file.filename+output_suffix)


    def make_binding_arguments(self, input_files: List[InputFile]) -> str:
        formatted_files_list = [f"-v {file.local_file_directory}:{file.docker_file_directory}" for file in input_files]
        return " ".join(formatted_files_list) + " "

    def format_input_files(self, executor: EXECUTION_TYPE, input_files: List[InputFile]) -> str:
        if executor == EXECUTION_TYPE.DOCKER:
            formatted_files_list = [f"{file.flag} {os.path.join(file.docker_file_directory, file.filename)}" for file in input_files] # uses local file directory
        else: # local executor
            formatted_files_list = [f"{file.flag} {os.path.join(file.local_file_directory, file.filename)}" for file in input_files] # uses docker file directory
        return " ".join(formatted_files_list)

    def extract_auxiliary_arguments(self, num_loci: int, num_cores: int) -> str:
        if num_loci is not None:
            num_loci_flag = "-e " + str(num_loci) + " "
        else:
            num_loci_flag = ""
        num_cores_flag = "-c " + str(num_cores) + " "
        return num_loci_flag + num_cores_flag


    def extract_run_command(self, start_point: STAGE, end_point: STAGE, file_arguments: FileArguments):
        not_implemented_messg = "MSMuTect does not implement this yet"
        if start_point == STAGE.BAMS:
            if end_point in [STAGE.MUTATIONS_FULL,  STAGE.MUTATIONS_EFFICIENT]:
                input_files = [file_arguments.tumor_file, file_arguments.normal_file, file_arguments.locus_file, file_arguments.output_file]
                run_flags = self.mutation_run_flags(end_point)
            elif end_point == STAGE.HISTOGRAMS or end_point == STAGE.ALLELES:
                input_files = [file_arguments.single_file, file_arguments.locus_file, file_arguments.output_file]
                if end_point == STAGE.ALLELES:
                    run_flags = "-A "
                else: # end_point == STAGE.HISTOGRAMS
                    run_flags = "-H "
        elif start_point == STAGE.HISTOGRAMS:
            run_flags = " --from_file "
            input_files = [file_arguments.precompiled_normal_hist_file, file_arguments.precompiled_tumor_hist_file,
                 file_arguments.output_file]
            if end_point in [STAGE.MUTATIONS_FULL, STAGE.MUTATIONS_EFFICIENT]:
                run_flags += self.mutation_run_flags(end_point)
            else:
                raise RuntimeError(not_implemented_messg)
        else:
            raise NotImplementedError(not_implemented_messg)
        return input_files, run_flags

    def mutation_run_flags(self, stage: STAGE):
        if stage == STAGE.MUTATIONS_FULL:
            return "-m -A -H "
        elif stage == STAGE.MUTATIONS_EFFICIENT:
            return "-m "
        else:
            raise RuntimeError("Is intended for mutation stage only")

if __name__ == "__main__":
    cmd_args = CMDArguments(EXECUTION_TYPE.DOCKER, STAGE.BAMS, STAGE.MUTATIONS_FULL)
    print(cmd_args.command)