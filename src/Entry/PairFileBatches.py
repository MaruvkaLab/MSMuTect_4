# cython: language_level=3
import os
from typing import List
from collections import namedtuple, defaultdict

from src.IndelCalling.Locus import Locus
from src.IndelCalling.AlleleSet import AlleleSet
from src.IndelCalling.Histogram import Histogram
from src.IndelCalling.CallAllelesFast import calculate_alleles
from src.IndelCalling import CallAllelesFast
from src.IndelCalling.CallMutations import call_mutations, is_possible_mutation
from src.IndelCalling.MutationCall import MutationCall

from src.GenomicUtils.ReadsFetcher import ReadsFetcher
from src.GenomicUtils.LocusFile import LociManager
from src.GenomicUtils.NoiseTable import get_noise_table

from src.Entry import BatchUtil
from src.Entry.FileBackedQueue import FileBackedQueue
from src.Entry.ResultsReaders import NamedColumns, to_int, str_to_bool
from src.Entry.InputHandler import exit_on

GERMLINE_VARIANT_COLUMN = "KNOWN_GERMLINE_VARIANT"

PairResults = namedtuple("PairResults", ['normal_alleles', 'tumor_alleles', 'decision'])



def format_mutation_call(decision: MutationCall):
    return f"{str(decision.normal_alleles.histogram.locus)}\t{str(decision.normal_alleles.histogram)}\t{str(decision.normal_alleles)}\t{str(decision.tumor_alleles.histogram)}\t{str(decision.tumor_alleles)}\t{str(decision)}"


def run_full_pair(normal: str, tumor: str, reference_genome_file: str, loci_file: str, batch_start: int,
                       batch_end: int, cores: int, flanking: int, required_reads: int, imprecise_mode: bool,
                  output_prefix: str) -> str:
    # returns path of output file
    loci_iterator = LociManager(loci_file, batch_start)
    noise_table = get_noise_table()
    results: List[FileBackedQueue] = BatchUtil.run_batch(partial_full_pair, [normal, tumor, reference_genome_file, flanking, noise_table, required_reads, imprecise_mode], loci_iterator,
                                  (batch_end - batch_start), cores, os.path.dirname(output_prefix))
    mutation_header = f"{Locus.header()}\t{Histogram.header(prefix='NORMAL_')}\t{AlleleSet.header(prefix='NORMAL_')}\t{Histogram.header(prefix='TUMOR_')}\t{AlleleSet.header(prefix='TUMOR_')}\t{MutationCall.header()}"
    output_file = output_prefix + ".full.mut"
    BatchUtil.write_queues_results(output_file, results, mutation_header)
    return output_file+".tsv"


def get_alleles(locus: Locus, reads_fetcher: ReadsFetcher, flanking: int, noise_table, required_reads: int,
                imprecise_mode: bool) -> AlleleSet:
    histogram = Histogram(locus, imprecise_mode)
    reads = reads_fetcher.get_reads(locus.chromosome, locus.start - flanking, locus.end + flanking)
    histogram.add_reads(reads)
    alleles = calculate_alleles(histogram, noise_table, required_read_support=required_reads)
    return alleles


def partial_full_pair(loci: List[Locus], normal: str, tumor: str, reference_genome_file: str, flanking: int,
                      noise_table, required_reads: int, imprecise_mode: bool,
                      results_dir: str) -> FileBackedQueue:
    calls = FileBackedQueue(out_file_dir=results_dir, max_memory=10**7)  # 10MB
    if len(loci) != 0:
        normal_fetcher = ReadsFetcher(normal, loci[0].chromosome, reference_genome_file)
        tumor_fetcher = ReadsFetcher(tumor, loci[0].chromosome, reference_genome_file)
        for locus in loci:
            normal_alleles = get_alleles(locus, normal_fetcher, flanking, noise_table, required_reads, imprecise_mode)
            tumor_alleles = get_alleles(locus, tumor_fetcher, flanking, noise_table, required_reads, imprecise_mode)
            calls.append(format_mutation_call(call_mutations(normal_alleles, tumor_alleles, noise_table)))
    calls.close()
    return calls


def run_mutations_pair(normal: str, tumor: str, reference_genome_file: str, loci_file: str, batch_start: int,
                       batch_end: int, cores: int, flanking: int, required_reads: int, imprecise_mode: bool,
                       output_prefix: str):
    # returns output file
    loci_iterator = LociManager(loci_file, batch_start)
    noise_table = get_noise_table()
    results: List[FileBackedQueue] = BatchUtil.run_batch(partial_mutations_pair, [normal, tumor, reference_genome_file, flanking, noise_table,
                                                                                  required_reads, imprecise_mode],
                                                     loci_iterator,
                                                     (batch_end - batch_start), cores, result_dir=os.path.dirname(output_prefix))
    mutation_header = f"{Locus.header()}\t{Histogram.header(prefix='NORMAL_')}\t{AlleleSet.header(prefix='NORMAL_')}\t{Histogram.header(prefix='TUMOR_')}\t{AlleleSet.header(prefix='TUMOR_')}\t{MutationCall.header()}"
    output_file = output_prefix + ".partial.mut"
    BatchUtil.write_queues_results(output_file, results, mutation_header)
    return output_file+".tsv"


def get_tumor_alleles(reads_fetcher: ReadsFetcher, locus: Locus, flanking: int, noise_table, required_reads=6, imprecise_mode=False) -> AlleleSet:
    histogram = Histogram(locus, imprecise_mode)
    reads = reads_fetcher.get_reads(locus.chromosome, locus.start - flanking, locus.end + flanking)
    histogram.add_reads(reads)
    current_alleles = calculate_alleles(histogram, noise_table, required_read_support=required_reads)
    return current_alleles


def partial_mutations_pair(loci: List[Locus], normal: str, tumor: str, reference_genome_file: str, flanking: int,
                           noise_table, required_reads: int, imprecise_mode: bool, results_dir: str) -> FileBackedQueue:
    calls = FileBackedQueue(out_file_dir=results_dir, max_memory=10**7) # 10MB
    if len(loci) != 0:
        normal_fetcher = ReadsFetcher(normal, loci[0].chromosome, reference_genome_file)
        tumor_fetcher = ReadsFetcher(tumor, loci[0].chromosome, reference_genome_file)
        for locus in loci:
            normal_alleles = get_alleles(locus, normal_fetcher, flanking, noise_table, required_reads, imprecise_mode)
            if is_possible_mutation(normal_alleles):
                tumor_alleles = get_alleles(locus, tumor_fetcher, flanking, noise_table, required_reads, imprecise_mode)
                calls.append(format_mutation_call(call_mutations(normal_alleles, tumor_alleles, noise_table)))
    calls.close()
    return calls


class HistogramTsvParser:
    """Builds Histogram objects from a histogram .tsv, addressing columns by name.

    Prefix-agnostic: works whether the motif/support columns are unprefixed
    (MOTIF_REPEATS_*, from single-file -H output) or prefixed (NORMAL_/TUMOR_,
    from a cut .full.mut.tsv). Robust to added/removed/reordered columns.
    """

    def __init__(self, header_line: str):
        self.columns = NamedColumns(header_line)
        self.motif_base = self.columns.resolve_base("MOTIF_REPEATS")
        self.support_base = self.columns.resolve_base("SUPPORTING_READS")
        self.has_germline_column = self.columns.has(GERMLINE_VARIANT_COLUMN)

    def parse(self, histogram_line: str) -> Histogram:
        fields = histogram_line.strip().split("\t")
        if len(fields) == 1:
            raise StopIteration("Iterated past end of file")
        c = self.columns
        locus = Locus(c.value(fields, "CHROMOSOME", str),
                      c.value(fields, "START", int),
                      c.value(fields, "END", int),
                      c.value(fields, "PATTERN", str),
                      c.value(fields, "REFERENCE_REPEATS", float),
                      c.value(fields, "REFERENCE_SEQUENCE", str),
                      known_germline_variant=c.value(fields, "KNOWN_GERMLINE_VARIANT", str_to_bool, default=False))
        motif_repeats = c.group(fields, self.motif_base, float) if self.motif_base else []
        motif_repeat_support = c.group(fields, self.support_base, to_int) if self.support_base else []
        repeat_dict = defaultdict(int)
        for repeat_length, support in zip(motif_repeats, motif_repeat_support):
            repeat_dict[int(repeat_length)] = int(support)
        histogram = Histogram(locus=locus)
        histogram.repeat_lengths = repeat_dict
        return histogram


def run_from_file(tumor_fp: str, normal_fp: str, batch_start: int, batch_end: int, required_reads: int,
                  output_prefix: str):
    noise_table = get_noise_table()
    results_dir = os.path.dirname(output_prefix)
    mutation_calls = FileBackedQueue(out_file_dir=results_dir, max_memory=int(1e7))
    tumor_file = open(tumor_fp, 'r')
    normal_file = open(normal_fp, 'r')
    tumor_parser = HistogramTsvParser(tumor_file.readline())  # consumes (and reads) the header line
    normal_parser = HistogramTsvParser(normal_file.readline())
    # The germline-variant annotation propagated to the output comes from the normal
    # histogram (see format_mutation_call), so require it to be present there.
    if not normal_parser.has_germline_column:
        exit_on(f"normal histogram file '{normal_fp}' has no {GERMLINE_VARIANT_COLUMN} column. "
                f"Regenerate it with a current version of MSMuTect (or carve it from a "
                f".full.mut.tsv that includes the {GERMLINE_VARIANT_COLUMN} field).")
    for i in range(batch_start - 1):
        tumor_file.readline()
        normal_file.readline()
    for i in range(batch_end - batch_start):
        try:
            tumor_histogram = tumor_parser.parse(tumor_file.readline())
        except StopIteration:
            break # finished consuming file. batch end could be malformed

        normal_histogram = normal_parser.parse(normal_file.readline())

        tumor_alleles = CallAllelesFast.calculate_alleles(tumor_histogram, noise_table, required_read_support=required_reads)
        normal_alleles = CallAllelesFast.calculate_alleles(normal_histogram, noise_table, required_read_support=required_reads)
        mutation_calls.append(format_mutation_call(call_mutations(normal_alleles, tumor_alleles, noise_table)))
        # exit()
    mutation_calls.close()
    mutation_header = f"{Locus.header()}\t{Histogram.header(prefix='NORMAL_')}\t{AlleleSet.header(prefix='NORMAL_')}\t{Histogram.header(prefix='TUMOR_')}\t{AlleleSet.header(prefix='TUMOR_')}\t{MutationCall.header()}"
    output_file = output_prefix + ".full.mut"
    BatchUtil.write_queues_results(output_file, [mutation_calls], mutation_header)
    tumor_file.close()
    normal_file.close()
    return output_file + ".tsv"


if __name__ == '__main__':
    # ../../ MSMuTect_0
    # .5 / msmutect.sh - -from_file - m - A - N
    # MSMuTect_normal_hist_to_run.normal.hist.tsv - T
    # MSMuTect_normal_hist_to_run.normal.hist.tsv - l
    # MSMuTect_normal_hist_to_run.normal.hist.tsv - O
    # res / croc
    #run_from_file(tumor_fp: str, normal_fp: str, batch_start: int, batch_end: int, required_reads: int, output_prefix: str)
    # tumor_fp = "/home/avraham/MaruvkaLab/msmutect_runs/gaia_from_file/MSMuTect_tumor_hist_to_run.tumor.hist.tsv"
    # normal_fp = "/home/avraham/MaruvkaLab/msmutect_runs/gaia_from_file/MSMuTect_normal_hist_to_run.normal.hist.tsv"
    # run_from_file(tumor_fp, normal_fp, 0, 10, 5, "/home/avraham/MaruvkaLab/msmutect_runs/gaia_from_file/res/croc")

    # tumor_fp = "/media/avraham/all_qs_data/tmp/0299c399-f301-444a-b7c4-57fa96ef5e2f.hist.tsv"
    # normal_fp = "/media/avraham/all_qs_data/tmp/0148fc66-05c6-489c-af47-26d12cbd898d.hist.tsv"

    tumor_fp = "/home/avraham/MaruvkaLab/msmutect_development/features/from_file_confirmation/tumor.hist.tsv"
    normal_fp = "/home/avraham/MaruvkaLab/msmutect_development/features/from_file_confirmation/normal.hist.tsv"

    run_from_file(tumor_fp, normal_fp, 0, 2000, 5,
                  "/home/avraham/MaruvkaLab/msmutect_development/features/from_file_confirmation/ff")
    # print(  f"{Locus.header()}\t{Histogram.header(prefix='NORMAL_')}\t{AlleleSet.header(prefix='NORMAL_')}\t{Histogram.header(prefix='TUMOR_')}\t{AlleleSet.header(prefix='TUMOR_')}\t{MutationCall.header()}")
# )
#     run_mutations_pair("/home/avraham/MaruvkaLab/Texas/texas_stad_run/tst/098698a0-3107-49e3-9226-d6d105f195a1.hist.tsv",
#                   "/home/avraham/MaruvkaLab/Texas/texas_stad_run/tst/009dcaf2-f6bb-415e-b088-6e852853b1a2.hist.tsv",
#                   44_347, 44_348, 5, True, "/home/avraham/MaruvkaLab/Texas/strict_msmutect/tmp")
    # run_from_file("/home/avraham/MaruvkaLab/Texas/texas_stad_run/tst/098698a0-3107-49e3-9226-d6d105f195a1.hist.tsv",
    #               "/home/avraham/MaruvkaLab/Texas/texas_stad_run/tst/009dcaf2-f6bb-415e-b088-6e852853b1a2.hist.tsv",
    #               16697, 16698, 5, True, "/home/avraham/MaruvkaLab/Texas/efficient_run/o16b")
    # run_from_file("/home/avraham/MaruvkaLab/Texas/texas_stad_run/tst/098698a0-3107-49e3-9226-d6d105f195a1.hist.tsv",
    #               "/home/avraham/MaruvkaLab/Texas/texas_stad_run/tst/009dcaf2-f6bb-415e-b088-6e852853b1a2.hist.tsv",
    #               16697, 16698, 5, True, "/home/avraham/MaruvkaLab/Texas/efficient_run/o16c")
    # run_from_file("/home/avraham/MaruvkaLab/Texas/texas_stad_run/tst/098698a0-3107-49e3-9226-d6d105f195a1.hist.tsv",
    #               "/home/avraham/MaruvkaLab/Texas/texas_stad_run/tst/009dcaf2-f6bb-415e-b088-6e852853b1a2.hist.tsv",
    #               16697, 16698, 5, True, "/home/avraham/MaruvkaLab/Texas/efficient_run/o16d")