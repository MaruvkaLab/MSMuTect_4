import unittest, subprocess, os, csv
from typing import Dict, List, Tuple

from tests.E2E_tests.load_config import bash_path, current_directory

# Carve a normal and a tumor histogram out of a .full.mut.tsv (selecting the locus
# columns plus the NORMAL_/TUMOR_ histogram columns by name), run
#     msmutect.sh -N <normal.hist> -T <tumor.hist> -m -A -O <prefix> --from_file
# and confirm the recomputed CALLs match the original .full.mut.tsv.

INPUT_DIR = os.path.join(current_directory(), "input_files")
OUTPUT_DIR = os.path.join(current_directory(), "output")

FULL_MUT_TSV = os.path.join(INPUT_DIR, "tst_cases_wfield.full.mut.tsv")
NORMAL_HIST = os.path.join(INPUT_DIR, "test_cases_normal.hist.tsv")
TUMOR_HIST = os.path.join(INPUT_DIR, "test_cases_tumor.hist.tsv")

# Histogram files are carved from the .full.mut.tsv by COLUMN NAME (not fixed
# positions), so the test is robust to schema changes (e.g. KNOWN_GERMLINE_VARIANT
# being added or Noisy Locus removed). Absent columns are simply skipped.
LOCUS_COLUMNS: List[str] = ["CHROMOSOME", "START", "END", "PATTERN",
                           "REFERENCE_SEQUENCE", "REFERENCE_REPEATS", "KNOWN_GERMLINE_VARIANT"]


def histogram_columns(prefix: str) -> List[str]:
    return ([f"{prefix}MOTIF_REPEATS_{i}" for i in range(1, 7)] +
            [f"{prefix}SUPPORTING_READS_{i}" for i in range(1, 7)])


Locus = Tuple[str, str, str]


class TestFromFilePreCalled(unittest.TestCase):

    def setUp(self):
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        self.ensure_histogram_files()

    def ensure_histogram_files(self):
        # Only (re)build the histograms when they are missing, mirroring the
        # `if [ ! -f ... ]` guard implied by the reference snippet.
        if os.path.exists(NORMAL_HIST) and os.path.exists(TUMOR_HIST):
            return
        self.cut_columns(FULL_MUT_TSV, LOCUS_COLUMNS + histogram_columns("NORMAL_"), NORMAL_HIST)
        self.cut_columns(FULL_MUT_TSV, LOCUS_COLUMNS + histogram_columns("TUMOR_"), TUMOR_HIST)

    @staticmethod
    def cut_columns(source: str, wanted_columns: List[str], destination: str):
        # select the wanted columns BY NAME from the source header; skip any not present
        with open(source, newline="") as src, open(destination, "w", newline="\n") as dst:
            reader = csv.reader(src, delimiter="\t")
            writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
            header = next(reader)
            index = {name: i for i, name in enumerate(header)}
            keep = [index[name] for name in wanted_columns if name in index]
            writer.writerow([header[i] for i in keep])
            for row in reader:
                writer.writerow([row[i] for i in keep])

    def test_full_from_file_calls_match_original(self):
        output_file = self.run_from_file("output_prefix_full")
        self.assert_calls_match_original(output_file)

    def test_efficient_from_file_calls_match_original(self):
        output_file = self.run_from_file("output_prefix_efficient")
        self.assert_calls_match_original(output_file)

    def run_from_file(self, output_prefix_name: str) -> str:
        output_prefix = os.path.join(OUTPUT_DIR, output_prefix_name)
        command = [
            bash_path(),
            "-N", NORMAL_HIST,
            "-T", TUMOR_HIST,
            "-m", "-A",
            "-O", output_prefix,
            "--from_file",
            "-f",  # allow the test to be re-run without manual cleanup
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode != 0:
            self.fail(f"from_file run {command} failed with return code {result.returncode}\n"
                      f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}")
        output_file = output_prefix + ".full.mut.tsv"
        self.assertTrue(os.path.exists(output_file), f"expected output file was not created: {output_file}")
        return output_file

    def assert_calls_match_original(self, output_file: str):
        original_calls = self.read_calls(FULL_MUT_TSV)
        produced_calls = self.read_calls(output_file)
        self.assertEqual(set(original_calls), set(produced_calls),
                         "loci in the from_file output differ from the original file")
        mismatches = {locus: (original_calls[locus], produced_calls[locus])
                      for locus in original_calls
                      if original_calls[locus] != produced_calls[locus]}
        self.assertEqual(mismatches, {},
                         f"{len(mismatches)} locus/loci got a different CALL than the original "
                         f"(showing up to 10): {dict(list(mismatches.items())[:10])}")

    @staticmethod
    def read_calls(path: str) -> Dict[Locus, str]:
        calls: Dict[Locus, str] = {}
        with open(path, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                locus: Locus = (row["CHROMOSOME"], row["START"], row["END"])
                calls[locus] = row["CALL"]
        return calls


if __name__ == "__main__":
    unittest.main()
