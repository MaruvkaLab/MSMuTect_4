#!/bin/bash

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PARENT_DIR="$( dirname "$SCRIPT_DIR" )"

PYTHONPATH="$PARENT_DIR"

# this hidden import system is insane, but other, cleaner approaches didn't work and it's not worth the time right now

pyinstaller --onefile --name msmutect \
--paths "$PARENT_DIR" \
--collect-all scipy \
--collect-all src.GenomicUtils \
--collect-submodules src.GenomicUtils \
--collect-submodules src.Entry \
--collect-submodules src.IndelCalling \
--hidden-import src.GenomicUtils.AlignmentFlags \
--hidden-import src.GenomicUtils.Mutation \
--hidden-import src.GenomicUtils.NoiseTable \
--hidden-import src.GenomicUtils.md_cigar_parser \
--hidden-import src.GenomicUtils.CigarOptions \
--hidden-import src.GenomicUtils.reference_locus_comparer \
--hidden-import src.IndelCalling.hist2vecs \
--hidden-import src.IndelCalling.AICs \
--hidden-import src.IndelCalling.CallAllelesFast \
--hidden-import src.IndelCalling.MutationCall \
--hidden-import src.Entry.FileBackedQueue \
--hidden-import src.Entry.InputHandler \
--hidden-import src.Entry.convert_tsv_to_vcf \
--hidden-import src.Entry.FormatUtil \
--hidden-import src.Entry.ResultsReaders \
--hidden-import src.Entry.BatchUtil \
--add-binary "$PARENT_DIR/src/GenomicUtils/*.so:src/GenomicUtils" \
--add-binary "$PARENT_DIR/src/IndelCalling/*.so:src/IndelCalling" \
--add-binary "$PARENT_DIR/src/Entry/*.so:src/Entry" \
"$PARENT_DIR/src/Entry/main.py"
