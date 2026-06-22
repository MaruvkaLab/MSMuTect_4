# cython: profile=True
# cython: language_level=3

import os
from typing import List, Optional
from pysam import AlignmentFile

from src.GenomicUtils.LocusFile import LociManager
from src.GenomicUtils.ReadsFetcher import ReadsFetcher
from src.GenomicUtils.NoiseTable import get_noise_table
from src.IndelCalling.Histogram import Histogram
from src.IndelCalling.AlleleSet import AlleleSet
from src.IndelCalling.Locus import Locus
from src.IndelCalling.CallAllelesFast import calculate_alleles
from src.Entry import BatchUtil
from src.Entry.FileBackedQueue import FileBackedQueue


def format_alleles(alleles: AlleleSet) -> str: # List[AlleleSet] not declared to avoid circular import
    return f"{alleles.histogram.locus}\t{str(alleles.histogram)}\t{str(alleles)}"


def run_single_allelic(BAM: str, reference_genome_file: str, loci_file: str, batch_start: int,
                       batch_end: int, cores: int, flanking: int, required_reads: int, imprecise_mode, output_prefix: str) -> None:
    loci_iterator = LociManager(loci_file, batch_start)
    noise_table = get_noise_table()
    results = BatchUtil.run_batch(partial_single_allelic, [BAM, reference_genome_file, flanking, noise_table, required_reads, imprecise_mode,],
                                                           loci_iterator,  (batch_end - batch_start), cores, os.path.dirname(output_prefix))
    header = f"{Locus.header()}\t{Histogram.header()}\t{AlleleSet.header()}"
    BatchUtil.write_queues_results(output_prefix + ".all", results, header)


def partial_single_allelic(loci: List[Locus], BAM: str, reference_genome_file: str, flanking: int, noise_table,
                           required_reads: int, imprecise_mode: bool, results_dir: str) -> FileBackedQueue:
    allelic_results = FileBackedQueue(out_file_dir=results_dir, max_memory=10**7) # 10MB
    if len(loci) != 0:
        reads_fetcher = ReadsFetcher(BAM, loci[0].chromosome, reference_genome_file)
        for locus in loci:
            current_histogram = Histogram(locus, imprecise_mode)
            reads = reads_fetcher.get_reads(locus.chromosome, locus.start - flanking, locus.end + flanking)
            current_histogram.add_reads(reads)
            current_alleles = calculate_alleles(current_histogram, noise_table, required_read_support=required_reads)
            allelic_results.append(format_alleles(current_alleles))
    allelic_results.close()
    return allelic_results


def run_single_histogram(BAM: str, reference_genome_file: str, loci_file: str, batch_start: int,
                         batch_end: int, cores: int, flanking: int, imprecise_mode: bool, output_prefix: str) -> None:
    loci_iterator = LociManager(loci_file, batch_start)
    results = BatchUtil.run_batch(partial_single_histogram, [BAM, reference_genome_file, flanking, imprecise_mode], loci_iterator,
                                  (batch_end - batch_start), cores, os.path.dirname(output_prefix))
    header = f"{Locus.header()}\t{Histogram.header()}"
    BatchUtil.write_queues_results(output_prefix + ".hist", results, header)


def format_histogram(histogram: Histogram) -> str:
    return f"{str(histogram.locus)}\t{str(histogram)}"


def partial_single_histogram(loci: List[Locus], BAM: str, reference_genome_file: str, flanking: int,
                             imprecise_mode: bool, results_dir: str) -> FileBackedQueue:
    histograms = FileBackedQueue(out_file_dir=results_dir, max_memory=10**7)  # 10MB
    if len(loci) != 0:
        reads_fetcher = ReadsFetcher(BAM, loci[0].chromosome, reference_genome_file)
        for locus in loci:
            current_histogram = Histogram(locus, imprecise_mode)
            reads = reads_fetcher.get_reads(locus.chromosome, locus.start - flanking, locus.end + flanking)
            current_histogram.add_reads(reads)
            histograms.append(format_histogram(current_histogram))
    histograms.close()
    return histograms



if __name__ == '__main__':
    # run_single_histogram(BAM: str, loci_file: str, batch_start: int,
    #                          batch_end: int, cores: int, flanking: int, output_prefix: str) -> None:
    # run_single_histogram("/home/avraham/MaruvkaLab/Texas/gdc/fc8376df-e7c3-42e9-b7b3-3413c3493874/02ca2f53-bd02-4e48-98be-aa09e69299ac_wgs_gdc_realn.bam",
    #                      "/home/avraham/MaruvkaLab/Texas/strict_msmutect/prob_locus.tmp",
    #                      0, 1, 1, 10, "/home/avraham/MaruvkaLab/Texas/strict_msmutect/tmp")
    #
    #

    run_single_histogram("/home/avraham/MaruvkaLab/msmutect_runs/data/other_bam/subset_N.bam",
                         "/home/avraham/MaruvkaLab/msmutect_runs/problems/seung_won_prob/phobos_hg19.validContig.sorted.txt/phobos_hg19.validContig.sorted.txt",
                         0, 10**12, 1, 10, "croc_trap")
    # run_single_histogram("/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/one_insertion_2_deletion.bam",
    #                  "/home/avraham/MaruvkaLab/MSMuTect_0.5/tests/sample_bams/real_locus.phobos",
    #                  0, 1, 1, 10, "/home/avraham/MaruvkaLab/bkup/garbage")
    # run_single_allelic(
    #     "/home/avraham/MaruvkaLab/msmutect_runs/problems/integer_wierd_rounding/prob_locus.bam",
    #     "/home/avraham/MaruvkaLab/msmutect_runs/problems/integer_wierd_rounding/prob_locus",
    #     0, 1, 1, 10, 6, True, "croc_tmp")
