import unittest, subprocess

class TestBatchUtil(unittest.TestCase):

    def test_docker_build(self):
        build_command = "docker buildx build --no-cache --output type=cacheonly ."
        self.run_build_command(build_command)

    def test_executable_build(self):
        build_command = "bash build_executable.sh"
        self.run_build_command(build_command)

    def run_build_command(self, build_command: str):
        results = subprocess.run(build_command, shell=True, capture_output=True)
        if results.returncode != 0:
            self.fail(f"{build_command} failed: {results.stderr}")



if __name__ == '__main__':
    unittest.main()
