import unittest, subprocess, os, glob

from tests.E2E_tests.CMDArguments import CMDArguments
from tests.E2E_tests.Enums import EXECUTION_TYPE, STAGE
from tests.E2E_tests.load_config import ground_truth_hist_file
from tests.testing_utils.high_level_results_parsing import count_mutations_in_file, entire_files_match


class TestBatchUtil(unittest.TestCase):
    # the full test suite would probably take around 10 hours to run. For a quick fix, just run the 100k tests, and test single file
    # HUNDRED_THOUSAND_LOCI_MUTATION_COUNT = 3377 # before new version of RR/LOH filtering
    HUNDRED_THOUSAND_LOCI_MUTATION_COUNT = 3095
    # HUNDRED_THOUSAND_LOCI_MUTATION_COUNT = 2925
    # FULL_RUN_MUTATION_COUNT = 832_982 # before new version of RR/LOH filtering
    FULL_RUN_MUTATION_COUNT = 831612

    # def test_from_file_100k_loci(self):
    #     arguments = CMDArguments(EXECUTION_TYPE.BASH, STAGE.HISTOGRAMS, STAGE.MUTATIONS_FULL, num_loci=100_000, num_cores=1)
    #     self.run_full_msmutect(arguments, self.HUNDRED_THOUSAND_LOCI_MUTATION_COUNT, remove=False)

    # def test_from_file_full(self):
    #     arguments = CMDArguments(EXECUTION_TYPE.BASH, STAGE.HISTOGRAMS, STAGE.MUTATIONS_FULL, num_cores=1)
    #     self.run_full_msmutect(arguments, self.FULL_RUN_MUTATION_COUNT, remove=False)

    # def test_single_file_full_executable(self):
    #     arguments = CMDArguments(EXECUTION_TYPE.EXECUTABLE, STAGE.BAMS, STAGE.HISTOGRAMS, num_cores=12)
    #     self.run_full_single_file(arguments, remove=False)

    # def test_full_loci_executable(self):
    #     arguments = CMDArguments(EXECUTION_TYPE.EXECUTABLE, STAGE.BAMS, STAGE.MUTATIONS_EFFICIENT, num_cores=16)
    #     self.run_full_msmutect(arguments, self.FULL_RUN_MUTATION_COUNT)

    # def test_full_loci_bash(self):
    #     arguments = CMDArguments(EXECUTION_TYPE.BASH, STAGE.BAMS, STAGE.MUTATIONS_EFFICIENT, num_cores=16)
    #     self.run_full_msmutect(arguments, self.FULL_RUN_MUTATION_COUNT, remove=False)

    # def test_100k_loci_docker_full(self):
    #     arguments = CMDArguments(EXECUTION_TYPE.DOCKER, STAGE.BAMS, STAGE.MUTATIONS_FULL, num_loci=100_000)
    #     self.run_full_msmutect(arguments, self.HUNDRED_THOUSAND_LOCI_MUTATION_COUNT)
    #
    # def test_100k_loci_executable_full(self):
    #     arguments = CMDArguments(EXECUTION_TYPE.EXECUTABLE, STAGE.BAMS, STAGE.MUTATIONS_FULL, num_loci=100_000)
    #     self.run_full_msmutect(arguments, self.HUNDRED_THOUSAND_LOCI_MUTATION_COUNT)
    #
    def test_100k_loci_bash_full(self):
        arguments = CMDArguments(EXECUTION_TYPE.BASH, STAGE.BAMS, STAGE.MUTATIONS_FULL, num_loci=100_000)
        self.run_full_msmutect(arguments, self.HUNDRED_THOUSAND_LOCI_MUTATION_COUNT)
    #
    # def test_100k_loci_docker_efficient(self):
    #     arguments = CMDArguments(EXECUTION_TYPE.DOCKER, STAGE.BAMS, STAGE.MUTATIONS_EFFICIENT, num_loci=100_000)
    #     self.run_full_msmutect(arguments, self.HUNDRED_THOUSAND_LOCI_MUTATION_COUNT)

    # #
    # def test_100k_loci_executable_efficient(self):
    #     arguments = CMDArguments(EXECUTION_TYPE.EXECUTABLE, STAGE.BAMS, STAGE.MUTATIONS_EFFICIENT, num_loci=100_000)
    #     self.run_full_msmutect(arguments, self.HUNDRED_THOUSAND_LOCI_MUTATION_COUNT)
    # #
    # def test_100k_loci_bash_efficient(self):
    #     arguments = CMDArguments(EXECUTION_TYPE.BASH, STAGE.BAMS, STAGE.MUTATIONS_EFFICIENT, num_loci=100_000)
    #     self.run_full_msmutect(arguments, self.HUNDRED_THOUSAND_LOCI_MUTATION_COUNT)

    def run_full_single_file(self, arguments: CMDArguments, remove: bool = True):
        self.run_command(arguments.command)
        self.assertTrue(entire_files_match(arguments.local_output_file, ground_truth_hist_file()))
        if remove:
            os.remove(arguments.local_output_file)

    def run_full_msmutect(self, arguments: CMDArguments, true_mutation_count: int, remove: bool = True):
        self.run_command(arguments.command)
        self.assert_correct_mutation_count(arguments.local_output_file, true_mutation_count)
        if remove:
            os.remove(arguments.local_output_file)

    def assert_correct_mutation_count(self, output_file_path: str, true_mutation_count: int):
        mutation_count = count_mutations_in_file(output_file_path)
        self.assertAlmostEqual(mutation_count, true_mutation_count, delta=int(true_mutation_count//10**2.5)) # bit random; 99.97% accurate

    def run_command(self, run_command: str):
        results = subprocess.run(run_command, shell=True, capture_output=True)
        if results.returncode != 0:
            self.fail( f"run_command {run_command} failed with return code {results.returncode}\nError: {results.stderr}")



if __name__ == '__main__':
    unittest.main()
