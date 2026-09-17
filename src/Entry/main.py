import argparse
import glob
import os
import random

from src.Entry.SingleFileBatches import run_single_allelic, run_single_histogram
from src.Entry.PairFileBatches import run_full_pair, run_mutations_pair, run_from_file
from src.Entry.InputHandler import create_parser, validate_input
from src.Entry.convert_tsv_to_vcf import convert_tsv_to_vcf
from src.Entry.format_additional_output_files import filter_mut_calls
from src.Entry.FileBackedQueue import RUN_ID_ENV_VAR, TMP_FILE_PREFIX


def count_lines(file: str):
    return sum(1 for _ in open(file, 'rb'))


def erase_tmp_files(output_prefix: str, run_id: str):
    # Temp files (see FileBackedQueue.get_unique_filename) are written to the
    # directory of the output prefix, named `<prefix><run_id>_...`. The run id is
    # shared by this process and every worker it forks, so scoping the glob to it
    # cleans up exactly this run's files and never a concurrent run's.
    output_dir = os.path.dirname(output_prefix) or os.getcwd()
    for tmp_file in glob.glob(os.path.join(output_dir, f"{TMP_FILE_PREFIX}{run_id}_*")):
        try:
            os.remove(tmp_file)
        except OSError:
            pass  # best-effort cleanup; don't mask the original error


def run_msmutect(args: argparse.Namespace):
    # Fix a per-run id before any worker forks so all of this run's temp files
    # (parent and workers) share it; workers inherit it via the environment.
    run_id = f"{os.getpid()}_{random.randint(0, 2 ** 31)}"
    os.environ[RUN_ID_ENV_VAR] = run_id
    try:
        _run_msmutect(args)
    except BaseException:
        # On any failure, remove temp files this run may have left behind.
        output_prefix = getattr(args, "output_prefix", None)
        if output_prefix:
            erase_tmp_files(output_prefix, run_id)
        raise


def _run_msmutect(args: argparse.Namespace):
    validate_input(args)  # will exit with error message if invalid combination of flags is given
    if args.from_file:
        if args.batch_end:
            batch_end = args.batch_end
        else:  # slight performance hit: ~ 1 sec / 2*10^6 loci
            batch_end = count_lines(args.tumor_file)
        mut_file = run_from_file(args.tumor_file, args.normal_file, args.batch_start, batch_end, args.read_level,
                                 args.output_prefix)
        # ----- fix -----
        filter_mut_calls(
            mut_file,
            call="M",
            call_column="CALL",
            sort_by="call"
        )
        if args.vcf:
            convert_tsv_to_vcf(mut_file, args.output_prefix + ".vcf")
    else:
        if args.batch_end:
            batch_end = args.batch_end
        else:  # slight performance hit: ~ 1 sec / 2*10^6 loci
            batch_end = count_lines(args.loci_file)
        if args.single_file:
            if args.allele or not args.histogram:
                run_single_allelic(args.single_file, args.reference_genome_file, args.loci_file, args.batch_start - 1,
                                   batch_end, args.cores, args.flanking, args.read_level, args.imprecise_mode, args.output_prefix)
            else:
                run_single_histogram(args.single_file, args.reference_genome_file, args.loci_file, args.batch_start - 1,
                                     batch_end, args.cores, args.flanking, args.imprecise_mode, args.output_prefix)

        else:
            if args.histogram and not args.mutation:
                run_single_histogram(args.normal_file, args.reference_genome_file, args.loci_file, args.batch_start - 1,
                                     batch_end, args.cores, args.flanking, args.imprecise_mode, args.output_prefix + ".normal")
                run_single_histogram(args.tumor_file, args.reference_genome_file, args.loci_file, args.batch_start - 1,
                                     batch_end, args.cores, args.flanking, args.imprecise_mode, args.output_prefix + ".tumor")
            elif args.allele and not args.mutation:
                run_single_allelic(args.normal_file, args.reference_genome_file, args.loci_file, args.batch_start - 1,
                                   batch_end, args.cores, args.flanking, args.read_level, args.imprecise_mode, args.output_prefix + ".normal")
                run_single_allelic(args.tumor_file, args.reference_genome_file, args.loci_file, args.batch_start - 1,
                                   batch_end, args.cores, args.flanking, args.read_level, args.imprecise_mode, args.output_prefix + ".tumor")
            else: # args.mutation=True
                if args.histogram or args.allele:
                    mut_file = run_full_pair(args.normal_file, args.tumor_file, args.reference_genome_file, args.loci_file, args.batch_start-1, batch_end,
                                  args.cores, args.flanking, args.read_level, args.imprecise_mode, args.output_prefix)

                else:  # args.mutation=True, just mutations
                    mut_file = run_mutations_pair(args.normal_file, args.tumor_file, args.reference_genome_file, args.loci_file, args.batch_start-1, batch_end,
                                    args.cores, args.flanking, args.read_level, args.imprecise_mode, args.output_prefix)
                # ---- fix -----
                filter_mut_calls(
                    mut_file,
                    call="M",
                    call_column="CALL",
                    sort_by="call"
                )
                if args.vcf:
                    convert_tsv_to_vcf(mut_file, args.output_prefix+".vcf")



if __name__ == "__main__":
    parser: argparse.ArgumentParser = create_parser()
    arguments = parser.parse_args()
    run_msmutect(arguments)
