import argparse, sys, os, pysam
from typing import List, Callable


def create_parser() -> argparse.ArgumentParser:
    # :return: creates parser with all command line arguments arguments
    MSMuTect_intro = "MSMuTect\n Version 4.1\n Authors: Yossi Maruvka, Avraham Kahan, and the Maruvka Lab at Technion. Please report bugs or issues to yosi.maruvka@bfe.technion.ac.il or k.avraham@technion.ac.il, or raise an issue on the github"
    parser = argparse.ArgumentParser(description=MSMuTect_intro)
    parser.add_argument("-T", "--tumor_file", help="Tumor BAM file")
    parser.add_argument("-N", "--normal_file", help="Non-tumor BAM file")
    parser.add_argument("-S", "--single_file", help="Analyze a single file for histogram and/or alleles")
    parser.add_argument("-l", "--loci_file", help="File of loci to be processed and included in the output", required=False)
    parser.add_argument("-O", "--output_prefix", help="prefix for all output files", required=True)
    parser.add_argument("-c", "--cores", help="Number of cores to run MSMuTect on", type=int, default=1)
    parser.add_argument("-b", "--batch_start", help="1-indexed number locus to begin analyzing at (Inclusive)", default=1, type=int)
    parser.add_argument("-e", "--batch_end", help="1-indexed number locus to stop analyzing at (Inclusive)", type=int)
    parser.add_argument("-H", "--histogram", help="Output a Histogram File", action='store_true')
    parser.add_argument("-A", "--allele", help="Output allele file", action='store_true')
    parser.add_argument("-m", "--mutation", help="Output mutation file", action='store_true')
    parser.add_argument("-F", "--flanking", help="Length of flanking on both sides of an accepted read", type=int, default=10)
    parser.add_argument("-r", "--read_level", help="Minimum number of reads to call allele", type=int, default=5)
    parser.add_argument("-f", "--force", help="Overwrite pre-existing files", action='store_true')
    parser.add_argument("--imprecise_mode", help="Avoids doing full local realignment against motif. Is much faster, but less precise, particularly for long motifs in impure loci", action='store_true')
    parser.add_argument("--reference_genome_file", help="If using CRAM files, you must provide a reference genome")
    parser.add_argument("--from_file", help="call alleles mutations on previously called indels (histograms). Only single core runs are supported", action='store_true')
    parser.add_argument("--vcf", help="Output VCF file in addition to TSV file. IMPORTANT NOTE: the location of the indel may be anywhere in the locus, but the VCF will always have it at the beginning of the locus. In addition, all loci with non-reference alleles are listed but only mutations are listed as 'PASS'", action='store_true')

    return parser


def exit_on(message: str, status: int = 1):
    # print message, and exit
    sys.stderr.write("ERROR: " + message + "\n")
    sys.exit(status)


def simple_index_check(bam: str):
    # Pick the right index extension and strip the correct file extension.
    if bam.endswith(".cram"):
        index_ext = ".crai"
        stem = bam[:-5]
    elif bam.endswith(".bam"):
        index_ext = ".bai"
        stem = bam[:-4]
    else:
        exit_on(f"Given file {bam} is not a .bam or .cram file")

    # Two accepted index naming conventions:
    #   samtools default:  file.bam.bai  / file.cram.crai
    #   alternative:       file.bai      / file.crai
    appended_index_path = bam + index_ext      # e.g. sample.bam.bai
    replaced_index_path = stem + index_ext     # e.g. sample.bai
    appended_file_exists = os.path.exists(appended_index_path)
    replaced_file_exists = os.path.exists(replaced_index_path)
    index_file_older_than_bam_message = f"Index file older than alignment file for {bam}. Index file must be younger than the alignment file. If you are sure the index file is correct, run 'touch [index_file]'"
    if appended_file_exists:
        if os.path.getmtime(appended_index_path) < os.path.getmtime(bam):
            exit_on(index_file_older_than_bam_message)
    elif replaced_file_exists:
        if os.path.getmtime(replaced_index_path) < os.path.getmtime(bam):
            exit_on(index_file_older_than_bam_message)
    else:
        exit_on(f"Given file {bam} is not sorted and/or indexed")


def validate_indexing(bam_files: List[str]) -> None:
    """ validates that given sequence files (BAM or CRAM) are indexed"""
    prefixes = ['', 'chr', 'Chr']
    for bam in bam_files:
        simple_index_check(bam)  # checks for the .bai / .crai file
        validated = False
        try:
            # Passing the path lets pysam auto-detect BAM vs CRAM from the file's
            # magic bytes (it opens CRAM in 'rc' mode automatically).
            current_handle = pysam.AlignmentFile(bam)
        except OSError:  # file was incomplete
            exit_on(f"{bam} or it's index file is incomplete")
        for prefix in prefixes:
            try:
                # CRAM file does not actually need reference file since it does not consume records
                _ = current_handle.fetch(f"{prefix}{1}", start=10_000, multiple_iterators=False)
                validated = True
                break  # verified
            except ValueError:  # different prefix, or not indexed
                continue
        if not validated:
            exit_on(f"Given file {bam} is not sorted and/or indexed, or contains an unusual prefix (not 'chr', 'Chr', or nothing)")



def validate_bams(arguments: argparse.Namespace):
    for f in [arguments.normal_file, arguments.tumor_file, arguments.single_file]:
        if not (f is None):
            if not os.path.exists(f):
                exit_on(f"Given sequencing file {f} does not exist")
            if not (f.endswith(".bam") or f.endswith(".cram")):
                exit_on(f"Given sequencing file {f} is not a BAM or CRAM file")
            if f.endswith(".cram"):
                if (arguments.reference_genome_file is None):
                    exit_on(f"When using CRAM files, you must provide a reference genome file")
                elif (not os.path.exists(arguments.reference_genome_file)):
                    exit_on(f"Given reference genome file {arguments.reference_genome_file} does not exist")
    if (bool(arguments.tumor_file) or bool(arguments.normal_file)) == bool(arguments.single_file): #  XOR
        exit_on("Provide Single file, or both Normal and Tumor file")
    elif bool(arguments.tumor_file) != bool(arguments.normal_file): #  XOR
        exit_on("Provide Single file, or both Normal and Tumor file")
    elif arguments.single_file:
        if not os.path.exists(arguments.single_file):
            exit_on("Provided single file path does not exist")
    else:
        if not os.path.exists(arguments.tumor_file) or not os.path.exists(arguments.normal_file):
            exit_on("Provided Normal or Tumor BAM path does not exist")
    validate_indexing([bam_file for bam_file in [arguments.tumor_file, arguments.normal_file, arguments.single_file] if bool(bam_file)])


def validate_output_files(arguments: argparse.Namespace):
    overwrite_files_mssg = "Files would be overwritten by this run. To force overwrite, use -f flag"
    if os.path.sep not in arguments.output_prefix:
        arguments.output_prefix = os.path.join(os.getcwd(), arguments.output_prefix)
    output_prefix = arguments.output_prefix
    if output_prefix.endswith(".tsv"):
        exit_on("Provided prefix ends with '.tsv'. The provided prefix is meant to be a a prefix, not the full output file path")
    if not os.path.exists(os.path.dirname(output_prefix)):
        exit_on("Output directory does not exist")
    if arguments.force:
        return
    elif arguments.single_file:
        if arguments.histogram and not arguments.allele:
            if os.path.exists(output_prefix + ".hist.tsv"):
                exit_on(overwrite_files_mssg)
        else:
            if os.path.exists(output_prefix + ".all.tsv"):
                exit_on(overwrite_files_mssg)
    else:  # pair file
        if arguments.mutation and arguments.vcf and os.path.exists(output_prefix+".vcf"):
            exit_on(overwrite_files_mssg)
        if (arguments.histogram or arguments.allele) and arguments.mutation and os.path.exists(output_prefix + ".full.mut.tsv"):
                    exit_on(overwrite_files_mssg)
        elif not arguments.histogram and not arguments.allele and os.path.exists(output_prefix + ".partial.mut.tsv"):
                exit_on(overwrite_files_mssg)
        elif arguments.allele and (os.path.exists(output_prefix + ".tumor.all.tsv") or
        os.path.exists(output_prefix + ".normal.all.tsv")):
            exit_on(overwrite_files_mssg)
        elif arguments.histogram and (os.path.exists(output_prefix + ".tumor.hist.tsv") or
                                      os.path.exists(output_prefix + ".normal.hist.tsv")):
            exit_on(overwrite_files_mssg)


def check_files(file_list: List[str], checking_function: Callable, template_error_message: str):
    for file_path in file_list:
        if not checking_function(file_path):
            exit_on(f"Problem with {file_path}: {template_error_message}")


def has_histogram_file_suffix(file_path: str) -> bool:
    return file_path.endswith(".hist.tsv")


def from_file_file_verification(arguments):
    if not (bool(arguments.tumor_file) and bool(arguments.normal_file)):
        exit_on("Provide both Normal and Tumor file for 'from_file' run")
    check_files([arguments.tumor_file, arguments.normal_file], os.path.exists, "file does not exist")
    check_files([arguments.tumor_file, arguments.normal_file], has_histogram_file_suffix, "file does not have proper suffix for a histogram file (.hist.tsv)")


def validate_input(arguments: argparse.Namespace):
    if arguments.from_file:
        if not arguments.mutation or arguments.single_file:
            exit_on("'from_file' option only supports calling mutations (since histograms were alreay called)")
        validate_output_files(arguments)
        from_file_file_verification(arguments)
    else:
        if not (arguments.loci_file or arguments.from_file): # there must be a source of loci; for from_file, this is the .hist.tsv files
            exit_on("Loci file must be provided")
        validate_bams(arguments)
        validate_output_files(arguments)
        if not os.path.exists(arguments.loci_file):
            exit_on("Provided loci file does not exist")
        if arguments.mutation and arguments.single_file:
            exit_on("Pair of files must be provided to call mutations")
        elif arguments.batch_start <= 0:
            exit_on("Batch Start must be equal to or greater than 1")
        elif arguments.cores <= 0:
            exit_on("Cores must be equal to or greater than 1")
        elif arguments.flanking < 0:
            exit_on("Flanking must be equal to or greater than 0")
        elif arguments.read_level < 1:
            exit_on("Minimum Read Level for calling alleles must be equal to or greater than 1")
        elif not os.path.exists(arguments.loci_file):
            exit_on("Loci file path does not exist")
        elif arguments.vcf and not arguments.mutation:
            exit_on("VCF file can only be generated for mutation calls")
