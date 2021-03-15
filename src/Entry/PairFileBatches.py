import numpy as np
from typing import List
from collections import namedtuple
from pysam import AlignmentFile

from src.IndelCalling.Locus import Locus
from src.IndelCalling.AlleleSet import AlleleSet
from src.IndelCalling.Histogram import Histogram
from src.IndelCalling.CallAlleles import calculate_alleles
from src.IndelCalling.CallMutations import call_mutations
from src.IndelCalling.FisherTest import Fisher
from src.GenomicUtils.ReadsFetcher import ReadsFetcher
from src.GenomicUtils.LocusFile import LociManager
from . import BatchUtil

PairResults = namedtuple("PairResults", ['normal_alleles', 'tumor_alleles', 'decision'])


def format_full_mutations(pairs: List[PairResults]) -> List[str]:
    output_lines = [f"{pair.normal_alleles.histogram.locus.chromosome}\t \
                    {pair.normal_alleles.histogram.locus.start}\t \
                    {pair.normal_alleles.histogram.locus.end}\t \
                    {pair.normal_alleles.histogram.locus.pattern}\t \
                    {pair.normal_alleles.histogram.locus.repeats}\t \
                    {pair.decision}\t \
                    {str(pair.normal_alleles.histogram)}\t \
                    {str(pair.tumor_alleles.histogram)}\t \
                    {str(pair.normal_alleles.log_likelihood)}\t \
                    {str(pair.tumor_alleles.log_likelihood)}\t \
                    {str(pair.normal_alleles.repeat_lengths)}\t \
                    {str(pair.tumor_alleles.repeat_lengths)}" for pair in pairs]
    return output_lines


def run_full_pair(normal: str, tumor: str, loci_file: str, batch_start: int,
                       batch_end: int, cores: int, flanking: int, output_prefix: str):
    loci_iterator = LociManager(loci_file, batch_start)
    noise_table = np.loadtxt(BatchUtil.get_noise_table_path(), delimiter=',')  # noise table
    results: List[PairResults] = BatchUtil.run_batch(partial_full_pair, [normal, tumor, flanking, noise_table], loci_iterator,
                                  (batch_end - batch_start) // 10_000, 10_000, cores)
    results += BatchUtil.run_batch(partial_full_pair, [normal, tumor, flanking, noise_table], loci_iterator, 1,
                                   (batch_end - batch_start) % 10_000, cores)
    allelic_header = "CHROMOSOME\tSTART\tEND\tPATTERN\tREPEATS\tHISTOGRAM\tLOG_LIKELIHOOD\tALLELES\tMUTATION_CALL"
    BatchUtil.write_results(output_prefix + ".normal.all", BatchUtil.format_alleles([pair.normal_alleles for pair in results]), allelic_header)
    BatchUtil.write_results(output_prefix + ".tumor.all", BatchUtil.format_alleles([pair.tumor_alleles for pair in results]), allelic_header)
    mutation_header = "CHROMOSOME\tSTART\tEND\tPATTERN\tREPEATS\tDECISION\tNORMAL_HISTOGRAM\tTUMOR_HISTOGRAM\tNORMAL_LOG_LIKELIHOOD\tTUMOR_LOG_LIKELIHOOD\tNORMAL_ALLELES\tTUMOR_ALLELES"
    BatchUtil.write_results(output_prefix + ".full.mut", format_full_mutations(results), mutation_header)


def get_allele_set(loci: List[Locus], BAM: str, flanking: int, noise_table) -> List[AlleleSet]:
    BAM_handle = AlignmentFile(BAM, "rb")
    reads_fetcher = ReadsFetcher(BAM_handle, loci[0].chromosome)
    allelic_results: List[AlleleSet] = []
    for locus in loci:
        current_histogram = Histogram(locus)
        reads = reads_fetcher.get_reads(locus.chromosome, locus.start - flanking, locus.end + flanking)
        current_histogram.add_reads(reads)
        current_alleles = calculate_alleles(current_histogram, noise_table)
        allelic_results.append(current_alleles)
    return allelic_results


def partial_full_pair(loci: List[Locus], normal: str, tumor: str, flanking: int, noise_table) -> List[PairResults]:
    if len(loci) == 0:
        return []
    normal_alleles: List[AlleleSet] = get_allele_set(loci, normal, flanking, noise_table)
    tumor_alleles: List[AlleleSet] = get_allele_set(loci, tumor, flanking, noise_table)
    fisher_calculator = Fisher()
    decisions = [call_mutations(normal_alleles[i], tumor_alleles[i], noise_table, fisher_calculator) for i in range(len(normal_alleles))]
    return [PairResults(normal_alleles=normal_alleles[i], tumor_alleles=tumor_alleles[i], decision = str(decisions[i]))
            for i in range(len(normal_alleles))]


def run_mutations_pair(tumor: str, normal: str, loci_file: str, batch_start: int,
                       batch_end: int, cores: int, flanking: int, output_prefix: str):
    pass


def run_single_allelic(BAM: str, loci_file: str, batch_start: int,
                       batch_end: int, cores: int, flanking: int, output_prefix: str) -> None:
    loci_iterator = LociManager(loci_file, batch_start)
    noise_table = np.loadtxt(BatchUtil.get_noise_table_path(), delimiter=',')  # noise table
    results = BatchUtil.run_batch(partial_single_allelic, [BAM, flanking, noise_table], loci_iterator,  (batch_end - batch_start)//10_000, 10_000, cores)
    results += BatchUtil.run_batch(partial_single_allelic, [BAM, flanking, noise_table], loci_iterator,  1, (batch_end - batch_start)%10_000, cores)
    header = "CHROMOSOME\tSTART\tEND\tPATTERN\tREPEATS\tHISTOGRAM\tLog_Likelihood\tALLELES"
    BatchUtil.write_results(output_prefix + ".all", results, header)


def partial_single_allelic(loci: List[Locus], BAM: str, flanking: int, noise_table) -> List[str]:
    BAM_handle = AlignmentFile(BAM, "rb")
    allelic_results: List[AlleleSet] = []
    if len(loci) == 0:
        return []
    reads_fetcher = ReadsFetcher(BAM_handle, loci[0].chromosome)

    return BatchUtil.format_alleles(allelic_results)