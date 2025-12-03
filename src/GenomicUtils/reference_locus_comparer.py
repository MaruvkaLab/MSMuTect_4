from dataclasses import dataclass
from typing import List, Dict, Tuple

from pysam.libcalignedsegment import AlignedSegment
from collections import namedtuple

from pysam.libcvcf import defaultdict

from src.GenomicUtils.CigarOptions import CIGAR_OPTIONS
from src.GenomicUtils.Mutation import Mutation
from src.GenomicUtils.md_cigar_parser import MD_OP, MD_TUPLE, read_md_ops_tuples


def relative_read_position(read_position: int, read_start: int, indel_bases: int):
    return read_position-read_start+indel_bases


def at_least_one_insertion(read: AlignedSegment) -> bool:
    at_least_one_insertion = False
    for cigar in read.cigartuples:
        cigar_op = cigar[0]
        if cigar_op == CIGAR_OPTIONS.INSERTION:
            return True
    return False

def extract_insertions(read: AlignedSegment, locus_start: int, locus_end: int):

    # returns all insertions WITHIN locus
    # locus end is inclusive
    read_start = read.reference_start + 1
    read_position = read_start
    indel_bases = 0
    if locus_end < read_position:
        return []

    if not at_least_one_insertion(read):
        return []

    mutations = []
    for cigar in read.cigartuples:
        cigar_op = cigar[0]

        cigar_op_len = cigar[1]
        if cigar_op in [CIGAR_OPTIONS.ALG_MATCH, CIGAR_OPTIONS.SEQ_MATCH, CIGAR_OPTIONS.SEQ_MISMATCH]:
            read_position += cigar_op_len

        elif cigar_op == CIGAR_OPTIONS.INSERTION:
            if locus_start<=read_position<=locus_end:
                relative_start = relative_read_position(read_position, read_start, indel_bases)
                new_indel = read.query_sequence[relative_start: relative_start + cigar_op_len]
                mutations.append(Mutation(new_indel, position=read_position, insertion=True))
            indel_bases+=cigar_op_len

        elif cigar_op == CIGAR_OPTIONS.DELETION:
            indel_bases -= cigar_op_len
            read_position += cigar_op_len

        if read_position > locus_end:
            break

    return mutations

def find_md_subsitutions_and_deletions(read: AlignedSegment, locus_start: int, locus_end: int,
                                       snp_padding) -> List[Mutation]:
    md_ops: List[MD_TUPLE] = read_md_ops_tuples(read)
    current_read_position = read.reference_start+1
    ret = []
    for op in md_ops:
        if op.op == MD_OP.MATCH:
            current_read_position+=op.length
        elif op.op == MD_OP.SUBSITUTION:
            if locus_start <= current_read_position <= locus_end or abs(locus_start - current_read_position) <= snp_padding or abs(
                    current_read_position - locus_end) <= snp_padding:
                ret.append(Mutation(op.seq, substitution=True, enters_or_exits_locus=False, position=current_read_position))
            current_read_position+=op.length # will always be 1
        elif op.op == MD_OP.DELETION:
            if locus_start <= current_read_position <= locus_end:
                # CHANGE: indel must start and end within locus
                if locus_end-current_read_position+1 >= op.length: # +1 for including the start and end
                    ret.append(Mutation(op.seq, position=current_read_position, deletion=True))
                else: # exits locus
                    ret.append(Mutation("IRRELEVANT", enters_or_exits_locus=True, deletion=True))
            elif current_read_position <= locus_start and locus_start<=(current_read_position+op.length): # enters locus
                ret.append(Mutation("IRRELEVANT", enters_or_exits_locus=True, deletion=True))
            current_read_position+=op.length

            # if locus_start <= current_read_position and locus_end-current_read_position+1 >= op.length: # +1 for including the start and end
            #     ret.append(Mutation(op.seq, position=current_read_position, deletion=True))
            # elif
            # current_read_position+=op.length
        else:
            raise RuntimeError("Failed to parse MD string")
        if current_read_position > (locus_end + snp_padding):
            break
    return ret


def extract_locus_mutations(read: AlignedSegment, locus_start: int, locus_end: int, snp_padding: int) -> List[Mutation]:

    # returns all indels WITHIN locus, and all snps within snp_padding bases of locus
    # locus end is inclusive
    insertions = extract_insertions(read, locus_start, locus_end)
    deletions_and_subsitutions = find_md_subsitutions_and_deletions(read, locus_start, locus_end, snp_padding)
    return insertions+deletions_and_subsitutions


def sort_str(a: str):
    return "".join(sorted(a))


def equivalent_strs_order_irrelevant(a: str, b: str):
    sorted_a = sort_str(a)
    sorted_b = sort_str(b)
    return sorted_a==sorted_b


@dataclass
class ReadContainers:
    # a read with an indel and a subsitution is an indel read
    full_match_reads: List[AlignedSegment]
    subsitution_reads: List[AlignedSegment]
    indel_reads: List[AlignedSegment]


def divide_reads_into_groups(reads: List[AlignedSegment]) -> ReadContainers:
    # find reads with subsitutions but no indels
    ret = ReadContainers([], [], [])
    for r in reads:
        all_cigar_ops = [cigar_op[0] for cigar_op in r.cigartuples]
        indel_read = False
        sub_read = False
        for op in all_cigar_ops:
            if op==CIGAR_OPTIONS.INSERTION or op==CIGAR_OPTIONS.DELETION:
                indel_read = True
                break
            elif op==CIGAR_OPTIONS.SEQ_MISMATCH:
                sub_read = True
        if indel_read:
            ret.indel_reads.append(r)
        elif sub_read:
            ret.subsitution_reads.append(r)
        else:
            ret.full_match_reads.append(r)
    return ret


def char_counts(pattern: str):
    char_count = defaultdict(int)
    for char in pattern:
        char_count[char]+=1
    return char_count

def equivalent_char_counts(a: str, b: str):
    a_char_count = char_counts(a)
    b_char_count = char_counts(b)
    for base in ["A", "C", "G", "T"]:
        if a_char_count[base]!=b_char_count[base]:
            return False
    return True

def twist(pattern: str, twist_slide: int):
    return pattern[twist_slide:] + pattern[:twist_slide]

def single_ms_indel_determination(pattern: str, indel: str) -> bool:
    if not equivalent_char_counts(pattern, indel):
        return False
    else:
        for i in range(len(pattern)):
            if twist(pattern, i) == indel:
                return True
        return False

def microsatellite_indel(indel: Mutation, ms_motif: str) -> int: # Tuple[int, bool]:
    # returns length of ms indel. if its not a valid ms indel, return 0
    # also returns whether there is an indel that eats into the proper locus
    if len(indel.seq)%len(ms_motif) != 0: # is not proper length
        return 0
    num_repeats = len(indel.seq)//len(ms_motif)
    if num_repeats==1:
        # could probably do something smarter to check (like a tree or whatever)... but not important enough to justify the effort
        if indel.deletion:
            factor = -1
        else: # insertion
            factor = 1
        return factor * int(single_ms_indel_determination(ms_motif, indel.seq))
    first_full_repeat_idx = indel.seq.find(ms_motif)
    if first_full_repeat_idx == -1:
        return 0
    else:
        reshuffled_locus = indel.seq[first_full_repeat_idx:]+indel.seq[:first_full_repeat_idx]
        if reshuffled_locus == ms_motif*num_repeats:
            length = num_repeats
        else:
            return 0

    if indel.insertion:
        return length
    else:
        return -length


if __name__ == '__main__':
    print(twist("ACTG", 0))
    print(twist("ACTG", 1))
    print(twist("ACTG", 2))
    print(twist("ACTG", 3))

    # a="AAAAAAAAAAAAAAAAAAAAAAAAACACACACACAAAAACGACGACGACGACAAAAAAAA"
    # print(char_count(a))
    # print(equivalent_strs_order_irrelevant("ACCT", "TCCDA"))
# @dataclass
# class ReadGroup:
#     seq: str
#     seq_char_count: str
#
# def group_equivalent_reads(reads: List[AlignedSegment]) -> Dict[str, List[AlignedSegment]]:
#     pass