# cython: language_level=3
import os
from typing import List
from collections import namedtuple, defaultdict
from pysam import AlignmentFile

from src.IndelCalling.Locus import Locus
from src.IndelCalling.AlleleSet import AlleleSet
from src.IndelCalling.Histogram import Histogram
from src.IndelCalling.CallAlleles import calculate_alleles
from src.IndelCalling import CallAllelesFast
from src.IndelCalling.CallMutations import call_mutations, is_possible_mutation
from src.IndelCalling.FisherTest import Fisher
from src.IndelCalling.MutationCall import MutationCall

from src.GenomicUtils.ReadsFetcher import ReadsFetcher
from src.GenomicUtils.LocusFile import LociManager
from src.GenomicUtils.NoiseTable import get_noise_table

from src.Entry import BatchUtil
from src.Entry.FileBackedQueue import FileBackedQueue

PairResults = namedtuple("PairResults", ['normal_alleles', 'tumor_alleles', 'decision'])


def format_mutation_call(decision: MutationCall):
    return f"{str(decision.normal_alleles.histogram.locus)}\t{str(decision.normal_alleles.histogram)}\t{str(decision.normal_alleles)}\t{str(decision.tumor_alleles.histogram)}\t{str(decision.tumor_alleles)}\t{str(decision)}"


def run_full_pair(normal: str, tumor: str, loci_file: str, batch_start: int,
                       batch_end: int, cores: int, flanking: int, required_reads: int, output_prefix: str) -> str:
    # returns path of output file
    loci_iterator = LociManager(loci_file, batch_start)
    noise_table = get_noise_table()
    results: List[FileBackedQueue] = BatchUtil.run_batch(partial_full_pair, [normal, tumor, flanking, noise_table, required_reads], loci_iterator,
                                  (batch_end - batch_start), cores, os.path.dirname(output_prefix))
    mutation_header = f"{Locus.header()}\t{Histogram.header(prefix='NORMAL_')}\t{AlleleSet.header(prefix='NORMAL_')}\t{Histogram.header(prefix='TUMOR_')}\t{AlleleSet.header(prefix='TUMOR_')}\t{MutationCall.header()}"
    output_file = output_prefix + ".full.mut"
    BatchUtil.write_queues_results(output_file, results, mutation_header)
    return output_file+".tsv"


def get_alleles(locus: Locus, reads_fetcher: ReadsFetcher, flanking: int, noise_table, required_reads: int) -> AlleleSet:
    histogram = Histogram(locus)
    reads = reads_fetcher.get_reads(locus.chromosome, locus.start - flanking, locus.end + flanking)
    histogram.add_reads(reads)
    alleles = calculate_alleles(histogram, noise_table, required_read_support=required_reads)
    return alleles


def partial_full_pair(loci: List[Locus], normal: str, tumor: str, flanking: int, noise_table, required_reads: int,
                      results_dir: str) -> FileBackedQueue:
    calls = FileBackedQueue(out_file_dir=results_dir, max_memory=10**7)  # 10MB
    if len(loci) != 0:
        normal_fetcher = ReadsFetcher(AlignmentFile(normal, "rb"), loci[0].chromosome)
        tumor_fetcher = ReadsFetcher(AlignmentFile(tumor, "rb"), loci[0].chromosome)
        fisher = Fisher()
        for locus in loci:
            normal_alleles = get_alleles(locus, normal_fetcher, flanking, noise_table, required_reads)
            tumor_alleles = get_alleles(locus, tumor_fetcher, flanking, noise_table, required_reads)
            calls.append(format_mutation_call(call_mutations(normal_alleles, tumor_alleles, noise_table, fisher)))
    calls.close()
    return calls


def run_mutations_pair(normal: str, tumor: str, loci_file: str, batch_start: int,
                       batch_end: int, cores: int, flanking: int, required_reads: int, output_prefix: str):
    # returns output file
    loci_iterator = LociManager(loci_file, batch_start)
    noise_table = get_noise_table()
    results: List[FileBackedQueue] = BatchUtil.run_batch(partial_mutations_pair, [normal, tumor, flanking, noise_table,
                                                                                  required_reads],
                                                     loci_iterator,
                                                     (batch_end - batch_start), cores, result_dir=os.path.dirname(output_prefix))
    mutation_header = f"{Locus.header()}\t{Histogram.header(prefix='NORMAL_')}\t{AlleleSet.header(prefix='NORMAL_')}\t{Histogram.header(prefix='TUMOR_')}\t{AlleleSet.header(prefix='TUMOR_')}\t{MutationCall.header()}"
    output_file = output_prefix + ".partial.mut"
    BatchUtil.write_queues_results(output_file, results, mutation_header)
    return output_file+".tsv"


def get_tumor_alleles(reads_fetcher: ReadsFetcher, locus: Locus, flanking: int, noise_table, required_reads=6) -> AlleleSet:
    histogram = Histogram(locus)
    reads = reads_fetcher.get_reads(locus.chromosome, locus.start - flanking, locus.end + flanking)
    histogram.add_reads(reads)
    current_alleles = calculate_alleles(histogram, noise_table, required_read_support=required_reads)
    return current_alleles


def partial_mutations_pair(loci: List[Locus], normal: str, tumor: str, flanking: int, noise_table, required_reads: int,
                           results_dir: str) -> FileBackedQueue:
    calls = FileBackedQueue(out_file_dir=results_dir, max_memory=10**7) # 10MB
    if len(loci) != 0:
        normal_fetcher = ReadsFetcher(AlignmentFile(normal, "rb"), loci[0].chromosome)
        tumor_fetcher = ReadsFetcher(AlignmentFile(tumor, "rb"), loci[0].chromosome)
        fisher = Fisher()
        for locus in loci:
            normal_alleles = get_alleles(locus, normal_fetcher, flanking, noise_table, required_reads)
            if is_possible_mutation(normal_alleles):
                tumor_alleles = get_alleles(locus, tumor_fetcher, flanking, noise_table, required_reads)
                calls.append(format_mutation_call(call_mutations(normal_alleles, tumor_alleles, noise_table, fisher)))
    calls.close()
    return calls


def construct_histogram_from_tsv(histogram_line: str) -> Histogram:
    histogram_line_cleaned = histogram_line.strip()
    broken_line = histogram_line_cleaned.split("\t")
    if len(broken_line) == 1:
        raise StopIteration("Iterated past end of file")
    dtypes = [str, int, int, str, str, float]
    chromosome, start, end, pattern, ref_seq, num_ref_repeats = (
    (dtypes[i](broken_line[i]) for i in range(len(dtypes))))
    motif_repeats = []
    motif_repeat_support = []
    for i in range(6, 12):
        if broken_line[i].strip() == "NA":
            break
        else:
            motif_repeats.append(float(broken_line[i]))
    for i in range(12, 18):
        if broken_line[i].strip() == "NA":
            break
        else:
            motif_repeat_support.append(int(broken_line[i]))
    locus = Locus(chromosome, start, end, pattern, num_ref_repeats, ref_seq)
    repeat_dict = defaultdict(int)
    for repeat_length, support in zip(motif_repeats, motif_repeat_support):
        repeat_dict[float(repeat_length)] = int(support)
    histogram = Histogram(locus=locus)
    histogram.repeat_lengths = repeat_dict
    return histogram


def run_from_file(tumor_fp: str, normal_fp: str, batch_start: int, batch_end: int, required_reads: int, output_prefix: str):
    noise_table = get_noise_table()
    fisher = Fisher()
    results_dir = os.path.dirname(output_prefix)
    mutation_calls = FileBackedQueue(out_file_dir=results_dir, max_memory=int(1e7))
    tumor_file = open(tumor_fp, 'r')
    normal_file = open(normal_fp, 'r')
    tumor_file.readline()  # burn header line
    normal_file.readline()  # burn header line
    for i in range(batch_start - 1):
        tumor_file.readline()
        normal_file.readline()
    for i in range(batch_end - batch_start):
        try:
            tumor_histogram = construct_histogram_from_tsv(tumor_file.readline())
        except StopIteration:
            break # finished consuming file. batch end could be malformed
        normal_histogram = construct_histogram_from_tsv(normal_file.readline())
        # tumor_alleles = calculate_alleles(tumor_histogram, noise_table,
        #                                                   required_read_support=required_reads)
        # normal_alleles = calculate_alleles(normal_histogram, noise_table,
        #                                                    required_read_support=required_reads)
        tumor_alleles = CallAllelesFast.calculate_alleles(tumor_histogram, noise_table, required_read_support=required_reads)
        normal_alleles = CallAllelesFast.calculate_alleles(normal_histogram, noise_table, required_read_support=required_reads)
        mutation_calls.append(format_mutation_call(call_mutations(normal_alleles, tumor_alleles, noise_table, fisher)))
        # exit()
    mutation_calls.close()
    mutation_header = f"{Locus.header()}\t{Histogram.header(prefix='NORMAL_')}\t{AlleleSet.header(prefix='NORMAL_')}\t{Histogram.header(prefix='TUMOR_')}\t{AlleleSet.header(prefix='TUMOR_')}\t{MutationCall.header()}"
    output_file = output_prefix + ".full.mut"
    BatchUtil.write_queues_results(output_file, [mutation_calls], mutation_header)
    tumor_file.close()
    normal_file.close()
    return output_file + ".tsv"


if __name__ == '__main__':
    print(  f"{Locus.header()}\t{Histogram.header(prefix='NORMAL_')}\t{AlleleSet.header(prefix='NORMAL_')}\t{Histogram.header(prefix='TUMOR_')}\t{AlleleSet.header(prefix='TUMOR_')}\t{MutationCall.header()}")
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