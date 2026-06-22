# cython: language_level=3
import random
from typing import Dict, List, Tuple
from collections import defaultdict

from pysam import AlignedSegment

from src.Entry.FormatUtil import format_list
from src.GenomicUtils.CigarOptions import CIGAR_OPTIONS
from src.GenomicUtils.Mutation import Mutation
from src.GenomicUtils.reference_locus_comparer import extract_locus_mutations, microsatellite_indel
from src.IndelCalling.Locus import Locus

#
# def static_vars(**kwargs): # thanks to stackoverflow John Kugelman
#     def decorate(func):
#         for k in kwargs:
#             setattr(func, k, kwargs[k])
#         return func
#     return decorate
#
# def return_minus_2():
#     return -1
#
#
# @static_vars(cache=defaultdict(float))
# def binom_cdf(n: int, k: int) -> float:
#
#     ret = binom_cdf.cache[]
#     if not in :
#         pass
#     else:
#         retur

class Histogram:
    def __init__(self, locus: Locus, imprecise_mode: bool = False):
        # imprecise mode means not to do a full local realignment against the actual motif
        # It's much faster, but less accurate, especially for long motifs and in impure loci
        self.imprecise_mode = imprecise_mode
        self.locus = locus
        self.repeat_lengths = defaultdict(int)  # key = repeat length; value = supporting reads
        self.noise_dict = defaultdict(int)
        self._noisiness = None

    def mutation_types(self, mutations: List[Mutation]) -> Tuple[bool, bool]:
        # returns whether at least one of the mutations is a snp (return value 1), or an indel (return value 2)
        has_snp = False
        has_indel = False
        for m in mutations:
            if m.substitution:
                has_snp = True
            elif m.insertion or m.deletion:
                has_indel = True
        return has_snp, has_indel

    def add_read_to_repeat_length_dict(self, read: AlignedSegment) -> int:
        mutations = extract_locus_mutations(read, self.locus.start, self.locus.end, len(self.locus.pattern))
        has_snp, has_indel = self.mutation_types(mutations)

        current_repeat_length = 0
        if has_snp:
            for m in mutations:
                if m.substitution:
                    self.noise_dict[m.position]+=1 # tracks how many snps there are at each location
            return # cannot support anything since it has a snp close to the locus
        elif has_indel:
            for mutation in mutations:
                if mutation.deletion or mutation.insertion:
                    if mutation.enters_or_exits_locus:
                        return
                    new_repeat_length = microsatellite_indel(mutation, self.locus.pattern)
                    current_repeat_length += new_repeat_length
            if current_repeat_length != 0:
                self.repeat_lengths[int(self.locus.repeats)+current_repeat_length]+=1
        else: # matches reference
            self.repeat_lengths[int(self.locus.repeats)]+=1

    def add_read_to_repeat_length_dict_efficient(self, read: AlignedSegment) -> float:
        # "old_mode": much faster, but also less accurate for longer motifs
        read_position = read.reference_start+1
        indel_bases = 0 # number of added/deleted bases in MS locus
        for cigar_op in read.cigartuples:

            if cigar_op[0] in [CIGAR_OPTIONS.ALG_MATCH, CIGAR_OPTIONS.SEQ_MATCH, CIGAR_OPTIONS.SEQ_MISMATCH]:
                read_position += cigar_op[1]
            elif cigar_op[0] == CIGAR_OPTIONS.INSERTION:
                if self.locus.start <= read_position <= self.locus.end:
                    indel_bases += cigar_op[1]
            elif cigar_op[0] == CIGAR_OPTIONS.DELETION:
                if read_position <= self.locus.end:
                    if read_position < self.locus.start:
                        deletion_length = max(cigar_op[1] + read_position - self.locus.start, 0)
                    else:
                        deletion_length = cigar_op[1]
                    indel_bases-=min(self.locus.end-read_position+1, deletion_length)
                read_position+=cigar_op[1]
        self.repeat_lengths[round(max(self.locus.repeats + indel_bases/len(self.locus.pattern), 0))]+=1 # so is never negative

    def add_reads(self, reads: List[AlignedSegment]) -> None:
        for read in reads:
            if self.imprecise_mode:
                self.add_read_to_repeat_length_dict_efficient(read)
            else: # default
                self.add_read_to_repeat_length_dict(read)


    def determine_if_locus_is_noisy(self):
        if len(self.noise_dict)==0: # no noise dict entries
            return False
        else:
            return max(self.noise_dict.values()) >= (max(5, sum(self.repeat_lengths.values())*0.3)) # should realy be replaced with binomial... problem is I'm not sure how to define the problem with mutation in founding cell of tumor

    def is_noisy(self):
        if self._noisiness is None:
            self._noisiness = self.determine_if_locus_is_noisy()
        return self._noisiness


    @property
    def rounded_repeat_lengths(self) -> defaultdict:
        # round all repeat lengths in histogram to nearest integer
        return self.repeat_lengths

    @staticmethod
    def header(prefix=''):
        return f"{prefix}MOTIF_REPEATS_1\t{prefix}MOTIF_REPEATS_2\t{prefix}MOTIF_REPEATS_3\t{prefix}MOTIF_REPEATS_4\t{prefix}MOTIF_REPEATS_5\t{prefix}MOTIF_REPEATS_6\t{prefix}SUPPORTING_READS_1\t{prefix}SUPPORTING_READS_2\t{prefix}SUPPORTING_READS_3\t{prefix}SUPPORTING_READS_4\t{prefix}SUPPORTING_READS_5\t{prefix}SUPPORTING_READS_6\tNoisy Locus"

    def prune_keys(self):
        for k in list(self.repeat_lengths.keys()): # list so dictionary size of keys don't change during pruning
            if self.repeat_lengths[k] == 0:
                del self.repeat_lengths[k]

    def __str__(self):
        self.prune_keys()
        sorted_repeats = sorted(self.repeat_lengths, key=self.repeat_lengths.get, reverse=True)
        ordered_repeats = [str(repeat) for repeat in sorted_repeats]
        ordered_support = [str(self.repeat_lengths[repeat]) for repeat in sorted_repeats]

        return format_list(ordered_repeats, 6) + "\t" + format_list(ordered_support, 6) + f"\t{int(self.is_noisy())}"

    def __eq__(self, other):
        for length in self.repeat_lengths:
            if not self.repeat_lengths[length] == other.repeat_lengths[length]:
                return False
        return len(self.repeat_lengths.keys()) == len(self.repeat_lengths.keys())

if __name__ == '__main__':
    print(Histogram.header())