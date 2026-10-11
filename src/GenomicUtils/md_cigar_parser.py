from dataclasses import dataclass
from typing import List, Tuple

from pysam import AlignedSegment


def split_cigar(cigar) -> List[str]:
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


class MD_OP:
    MATCH=1
    SUBSITUTION=2
    DELETION=3
    INSERTION=4


@dataclass
class MD_TUPLE:
    op: int
    length: int
    seq: str

def split_cigar_into_ops(cigar: str) -> List[Tuple[MD_OP, int]]:
    cigar_tuple_strs = split_cigar(cigar)
    cigar_translations = {"M": MD_OP.MATCH, "I": MD_OP.INSERTION, "X": MD_OP.SUBSITUTION, "D": MD_OP.DELETION}
    return [(cigar_translations[cig[-1]], int(cig[:-1])) for cig in cigar_tuple_strs]



def split_md_string_into_tokens(md_str: str) -> List[str]:
    # splits md string into tokens (ints, ops, letters)
    ret = []
    current_num = []
    for char in md_str:
        if str.isnumeric(char):
            current_num.append(char)
        else:
            if len(current_num)> 0:
                ret.append("".join(current_num))
                current_num=[]
            ret.append(char)
    if len(current_num)>0:
        ret.append("".join(current_num))
    return ret

def split_MD_into_ops(md: str) -> List[MD_TUPLE]:
    is_in_deletion=False
    current_deletion_bases = []
    ret = []
    md_tokens = split_md_string_into_tokens(md)
    for token in md_tokens:
        if token=='^':
            is_in_deletion = True
        elif str.isalpha(token):
            if is_in_deletion:
                current_deletion_bases.append(token)
            else: # otherwise is just part of deletion
                ret.append(MD_TUPLE(MD_OP.SUBSITUTION, 1, "N"))
        else: # char == num
            if is_in_deletion:
                ret.append(MD_TUPLE(MD_OP.DELETION, len(current_deletion_bases), "".join(current_deletion_bases)))
                current_deletion_bases=[]
                is_in_deletion=False
            ret.append(MD_TUPLE(MD_OP.MATCH, int(token), "N"))
    if is_in_deletion:
        ret.append(MD_TUPLE(MD_OP.DELETION, len(current_deletion_bases), "".join(current_deletion_bases)))
    return ret


def split_MD_into_tuples(md: str) -> List[Tuple[int, int]]:
    is_in_deletion=False
    current_deletion_length = 0
    ret = []
    md_tokens = split_md_string_into_tokens(md)
    for token in md_tokens:
        if token=='^':
            is_in_deletion = True
        elif str.isalpha(token):
            if is_in_deletion:
                current_deletion_length+=1
            else: # otherwise is just part of deletion
                ret.append((MD_OP.SUBSITUTION, 1))
        else: # char == num
            if is_in_deletion:
                ret.append((MD_OP.DELETION, current_deletion_length))
                current_deletion_length=[]
                is_in_deletion=False
            ret.append((MD_OP.MATCH, int(token)))
    if is_in_deletion:
        ret.append((MD_OP.DELETION, current_deletion_length))
    return ret

def extract_md_str(read: AlignedSegment) -> str:
    tags = read.get_tags()
    for t in tags:
        if t[0]=="MD":
            return t[1]
    raise RuntimeError("No MD tag found. MSMuTect needs MD tag to run. Earlier versions of MSMuTect (4.0 and earlier) can run without it")

def read_md_ops_tuples(read: AlignedSegment):
    md_str = extract_md_str(read)
    md_ops = split_MD_into_ops(md_str)
    return md_ops

# def unify_md_and_cigar_strs(cigar: str, md: str) -> str:
#     ret_ops = []
#     op_translation = {MD_OP.MATCH: "M", MD_OP.SUBSITUTION: "X", MD_OP.DELETION: "D", MD_OP.INSERTION: "I"}
#     md_ops = split_MD_into_ops(md)
#     cigar_ops = split_cigar_into_ops(cigar)
#     cigar_op_idx = 0
#     for md_idx in range(len(md_ops)):
#         current_md = md_ops[md_idx]
#         current_cigar = cigar_ops[cigar_op_idx]
#         current_cigar_length = current_cigar[1]
#         if current_md[0] == MD_OP.MATCH and current_cigar[0] == MD_OP.MATCH and current_cigar[1] < current_md[1]:
#             last_cigar_length = current_cigar_length
#             ret_ops.append(str(current_cigar_length)+op_translation[current_cigar_length])
#             cigar_op_idx+=1
#             current_cigar = cigar_ops[cigar_op_idx]
#             current_cigar_length = current_cigar[1]
#             assert current_cigar[0] == MD_OP.INSERTION
#             ret_ops.append(str(current_cigar_length) + op_translation[current_cigar[1]])
#             cigar_op_idx += 1
#             current_cigar = cigar_ops[cigar_op_idx]
#             current_cigar_length = current_cigar[1]-last_cigar_length # to adjust M
#         elif current_cigar:
#             ret_ops.append()
#     return "".join(ret_ops)






if __name__ == '__main__':
    print(split_cigar("78M54D3I1222I"))
