import unittest, subprocess
from pathlib import Path

# Marker files that only exist at the top-level MSMuTect directory. The build
# commands need these, so we locate the root by their presence.
_ROOT_MARKERS = ("Dockerfile", "build_executable.sh")


def _find_repo_root() -> Path:
    """Walk up from this file until we find the MSMuTect root.

    Anchoring on marker files (rather than a fixed number of parent hops)
    keeps this working even if the test file is moved to a different depth.
    """
    for directory in (Path(__file__).resolve(), *Path(__file__).resolve().parents):
        if all((directory / marker).is_file() for marker in _ROOT_MARKERS):
            return directory
    raise RuntimeError(
        f"Could not locate MSMuTect root (searched upward from {__file__} "
        f"for {_ROOT_MARKERS})"
    )


# Top-level MSMuTect directory, resolved relative to this file so the build
# commands work no matter which directory the tests are launched from.
REPO_ROOT = _find_repo_root()

class TestBatchUtil(unittest.TestCase):

    def setUp(self):
        pass

    def test_docker_build(self):
        build_command = "docker buildx build --no-cache --output type=cacheonly ."
        self.run_build_command(build_command)

    # def test_executable_build(self):
    #     build_command = "bash build_executable.sh"
    #     self.run_build_command(build_command)
    #
    def run_build_command(self, build_command: str):
        results = subprocess.run(build_command, shell=True, capture_output=True, cwd=REPO_ROOT)
        if results.returncode != 0:
            self.fail(f"{build_command} failed: {results.stderr}")



if __name__ == '__main__':
    unittest.main()
