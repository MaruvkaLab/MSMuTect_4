import re, os, shutil
from typing import List, Tuple

from dataclasses import dataclass

from tests.testing_utils.sample_sequences import simple_seq
from tests.testing_utils.self_contained_utils import header_only_sam, sample_bams_path


@dataclass
class FakeRead:
    read_start: int
    cigar_str: str
    sequence: str = None
    insertions_snps: List[str] = None

def split_cigar(cigar: str):
    type_idxs = []
    for i, c in enumerate(cigar):
        if c.isalpha():
            type_idxs.append(i)
    last_idx = 0
    ret = []
    for t in type_idxs:
        ret.append(cigar[last_idx:t+1])
        last_idx=t+1
    return ret

def create_MD_string(read: FakeRead, modifications: List[str]) -> str:
    """
    Creates MD string based on FakeRead and substitutions
    :param read: has cigar info, sequence, etc.
    :param modifications: Includes deletions!
    :return: the md string
    """
    cigar_split = split_cigar(read.cigar_str)
    ops = []
    sub_idx = 0
    current_match = 0
    for cig in cigar_split:
        if cig[-1]=="M":
            current_match+=int(cig[:-1])
            continue
        else:
            if current_match!=0:
                ops.append(current_match)
            current_match=0

        if cig[-1]=="X":
            ops.append(modifications[sub_idx])
            sub_idx+=1
        elif cig[-1] == "D":
            ops.append(f"^{modifications[sub_idx]}")
            sub_idx+=1

    if current_match!=0:
        ops.append(str(current_match))
    ret = [ops[0]]
    for i in range(1,len(ops)): # combine matches in case of MATCH-INSERTION-MATCH. Not the cleanest code
        if type(ret[-1]) == int and type(ops[i])==int:
            ret[-1]+=ops[i]
        else:
            ret.append(ops[i])
    ret = [str(r) for r in ret]
    return "MD:Z:"+"".join(ret)

def write_seq(read_start: int, cigar_str: str, insertions_snps: List[str], base_sequence: str, base_sequence_position=9975) -> Tuple[str, List[str]]:
    """
    returns properly formatted DNA sequence with given mutations (substitutions, indels)
    :param read_start: read start position (absolute)
    :param cigar_str: CIGAR string of read
    :param insertions_snps: insertions+substitutions sequences
    :param base_sequence: the base sequence to use
    :param base_sequence_position: the absolute posititon of the base sequence
    :return: (sequence, substitutions)
    """
    if read_start < base_sequence_position:
        raise RuntimeError
    op_lens = [int(op_len) for op_len in re.split("[MXID]", cigar_str)[:-1]]
    ops = [char for char in cigar_str if char in ["M", "D", "I", "X"]]
    current_pos = read_start-base_sequence_position
    segments = []
    if insertions_snps is None:
        insertions_snps = []
    modified_substitutions = insertions_snps.copy()
    sub_idx=0
    for op, op_len  in zip(ops, op_lens):
        if op == 'M':
            segments.append(base_sequence[current_pos:current_pos+op_len])
            current_pos+=op_len
        elif op == 'X':
            if insertions_snps[sub_idx] == "N": # correct if user doesn't bother setting
                if base_sequence[current_pos]=="G":
                    modified_substitutions[sub_idx] = "T"
                else:
                    modified_substitutions[sub_idx] = "G"
            assert op_len==1
            assert insertions_snps[sub_idx] != base_sequence[current_pos]
            segments.append(insertions_snps[sub_idx])
            sub_idx+=1
            current_pos+=1
        elif op=='D':
            modified_substitutions.append(base_sequence[current_pos:current_pos + op_len])
            current_pos+=op_len
        else: # I
            segments.append(insertions_snps[sub_idx])
            sub_idx+=1
    return "".join(segments), modified_substitutions


def add_to_end_of_file(fp: str, lines: str):
    with open(fp, 'a') as sam:
        sam.write(lines)

def create_readline(fake_reads: List[FakeRead], base_seq=None, base_seq_position: int = None, read_length: int = 101):
    if base_seq is None:
        base_seq = simple_seq()
    if base_seq_position is None:
        base_seq_position = 9975
    if len(fake_reads) == 0:
        return
    seq = 'A' * read_length

    # based on real TCGA case, but fields obscured and changed
    fields = ["FAKE",  # name
              '2', # bitwise flags. 2 indicates aligned properly
              '1', # chromosome
              '10051',  # start
              '4',  # mapping quality
              '101M',  # cigar
              '=', # reference name of next read
              '44444', # position of next read
              '108', # template length
              seq,
              seq,
              'X0:i:100',
              'MD:Z:79A22',
              'RG:Z:0.3',
              'XG:i:0',
              'AM:i:0',
              'NM:i:1',
              'SM:i:0',
              'XM:i:1',
              'XO:i:0',
              'MQ:i:0',
              'OQ:Z:C@SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS?A=A;ABB((,:9<??BBB@D9<A?B99AB',
              "XT:A:R"]


    all_new_reads = []
    for fr in fake_reads:
        new_read = fields.copy()
        new_read[3] = str(fr.read_start)
        new_read[5] = fr.cigar_str
        if fr.sequence is not None:
            new_read[9] = fr.sequence
            new_read[10] = fr.sequence
        else:
            new_read[9], updated_subsitutions = write_seq(fr.read_start, fr.cigar_str, fr.insertions_snps, base_seq, base_seq_position)
            new_read[10], _ = write_seq(fr.read_start, fr.cigar_str, fr.insertions_snps, base_seq, base_seq_position)

        new_read[12] = create_MD_string(fr, updated_subsitutions)

        all_new_reads.append(new_read)
    return "\n".join(["\t".join(a) for a in all_new_reads])

def create_new_bam(new_name: str, fake_reads: List[FakeRead], base_seq=None, base_seq_pos=None, read_length: int = 101):
    header_only_file = header_only_sam()
    new_filename = os.path.join(sample_bams_path(), new_name + ".sam")
    shutil.copyfile(header_only_file, new_filename)
    readlines = create_readline(fake_reads, base_seq, base_seq_pos, read_length)
    add_to_end_of_file(new_filename, readlines)
    bam_filename = new_filename[:-4]+".bam"
    os.system(f"samtools view -b -h {new_filename} > {bam_filename}")
    os.system(f"samtools index {bam_filename}")