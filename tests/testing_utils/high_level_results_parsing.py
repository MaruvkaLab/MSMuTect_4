import csv, os


def count_mutations_in_file(call_filepath: str) -> int:
    # returns how many mutations there are in a .full.mut.tsv or .partial.mut.tsv
    count = 0
    with open(call_filepath, mode="r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("CALL") == "M":
                count += 1
    return count

def entire_files_match(path1: str, path2: str, chunk_size: int = 8192) -> bool:
    """
    Return True if the two files match byte-for-byte, False otherwise.
    Reads in chunks to avoid loading entire files into memory.
    """
    # Quick check: if sizes differ, files can't be identical
    if os.path.getsize(path1) != os.path.getsize(path2):
        return False

    with open(path1, "rb") as f1, open(path2, "rb") as f2:
        while True:
            b1 = f1.read(chunk_size)
            b2 = f2.read(chunk_size)

            if not b1 and not b2:
                return True
    return True

if __name__ == "__main__":
    print(count_mutations_in_file("/home/avraham/MaruvkaLab/MSMuTect/tests/E2E_tests/output/test_output.partial.mut.tsv"))