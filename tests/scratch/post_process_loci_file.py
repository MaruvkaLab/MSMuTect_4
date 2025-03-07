import Levenshtein
# from difflib import SequenceMatcher


def reverse_complement(seq: str) -> str:
    complement_dict = ["T", "G", "C", "A"] # reverse of A,C,G,T
    seq_length = len(seq)
    ret_seq = ["N" for i in range(seq_length)]
    for i in range(len(seq)):
        current_base = seq[i]
        complement_idx = (current_base=='C')+2*(current_base=="G")+3*(current_base=="T")
        ret_seq[seq_length-i-1] = complement_dict[complement_idx]
    return "".join(ret_seq)

def extrapolate_pattern(pattern, dst_len):
    num_multiples = dst_len//len(pattern)
    partial_repeat = dst_len%len(pattern)
    ret = pattern*num_multiples+pattern[:partial_repeat]
    return ret


def similarity(a: str, b: str):
    # return SequenceMatcher(None, a, b).ratio()
    return Levenshtein.ratio(a, b)

def select_best_pattern(sequence: str, pattern: str) -> str:
    if pattern=="AGC":
        croc=1
    forward_pattern_seq = extrapolate_pattern(pattern , len(sequence))
    reverse_pattern_seq = reverse_complement(forward_pattern_seq)
    if sequence == reverse_pattern_seq:
        return reverse_complement(pattern)
    elif sequence == forward_pattern_seq:
        return pattern
    else:
        forward_similarity = similarity(sequence, forward_pattern_seq)
        reverse_similarity = similarity(sequence, reverse_pattern_seq)
        if max(forward_similarity, reverse_similarity) < 0.8 and 0.5 > abs(forward_similarity-reverse_similarity):
            # print(f"Problem: {pattern}, {sequence}, {forward_similarity>reverse_similarity}")
            croc=1
        if forward_similarity > reverse_similarity:
            return pattern
        else:
            return reverse_complement(pattern)

def create_new_loci_file(input_fp: str, output_fp: str):
    with open(input_fp, 'r') as input_loci:
        with open(output_fp, 'w+') as output_loci:
            current_line = input_loci.readline()
            split_line = current_line.split("\t")
            while current_line != "":
                split_line[12] = select_best_pattern(split_line[13], split_line[12])
                output_loci.write("\t".join(split_line))
                current_line = input_loci.readline()
                split_line = current_line.split("\t")


if __name__ == '__main__':
    # print(reverse_complement("ACGGGT"))
    # print(extrapolate_pattern("ACTG", 9))
    create_new_loci_file("/home/avraham/MaruvkaLab/Texas/SNVs/slimeball/data/hg38_1to15_all_perf_sorted.phobos",
                         "/home/avraham/MaruvkaLab/msmutect_runs/hg38_all_rc_corrected")