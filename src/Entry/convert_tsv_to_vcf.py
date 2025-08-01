import shutil
from typing import List
from src.Entry.ResultsReaders import ResultsReaderMutationFile, LocusMutationCall


def generate_ref_str_alternatives(reference: str, microsatellite: str, alternative_num_repeats: List[float]):
    ms_length = len(microsatellite)
    infinitely_long_ms = 1000*microsatellite
    alt_strs = []
    locus_length = len(reference)
    for alt in alternative_num_repeats:
        num_bases = int(ms_length*alt)
        num_added_bases = num_bases-locus_length
        if alt == 0:
            alt_strs.append("")
        else:
            if num_added_bases < 0: # deletion in alt
                alt_strs.append(reference[-1*num_added_bases:])
            elif num_added_bases > 0: # insertion in alt
                alt_strs.append(infinitely_long_ms[:num_bases]+reference)
    return alt_strs


def formulate_vcf_line(current_line: LocusMutationCall) -> str:
    # ref_repeats = current_line.num_ref_repeats
    alt_strs_normal = set(generate_ref_str_alternatives(current_line.ref_seq, current_line.pattern, current_line.normal_motif_repeats))
    alt_strs_tumor = set(generate_ref_str_alternatives(current_line.ref_seq, current_line.pattern, current_line.tumor_motif_repeats))
    all_alt_strs = alt_strs_tumor | alt_strs_normal

    output_line = [
        current_line.chromosome,
        str(current_line.start),
        ".", # id
        current_line.ref_seq,
        ",".join(all_alt_strs),
        "." # quality,
        "PASS", # FILTER
        current_line.mutation_call, # info
    ]
    if current_line.mutation_call.strip() == "M":
        filter = "PASS"
    else:
        filter = current_line.mutation_call

    output_line = [
        current_line.chromosome,
        str(current_line.start),
        ".",  # id
        current_line.ref_seq,
        ",".join(all_alt_strs),
        ".",  # quality,
        filter,  # FILTER
        "."  # info
    ]
    return "\t".join(output_line)


def create_vcf_lines(input_tsv_fp: str) -> str:
    output_lines = []
    results_reader = ResultsReaderMutationFile(input_tsv_fp)
    current_line = next(results_reader, None)
    while current_line is not None:
        num_normal_alleles = len(current_line.normal_motif_repeats)
        num_tumor_alleles = len(current_line.tumor_motif_repeats)
        if num_normal_alleles > 1 or num_tumor_alleles > 1:
            new_line = formulate_vcf_line(current_line)
            output_lines.append(new_line)
        current_line = next(results_reader, None)

    return "\n".join(output_lines)


def convert_tsv_to_vcf(input_tsv_fp: str, output_vcf_fp: str):
    vcf_lines = create_vcf_lines(input_tsv_fp)
    header = """\
##fileformat=VCFv4.0
##<ID=ID,Number=,Type=,Description="NOTE: FOR ALL IMPURE LOCI, THE LOCATION OF THE MUTATION MAY BE INCORRECT BY UP TO THE LENGTH OF THE MICROSATELLITE LOCUS. THE EXACT REF AND ALT SEQUENCE MAY BE ERRANT AS WELL">
##fileformat=VCFv4.2
##source=MSMuTect
##phasing=partial
##contig=<ID=1>
##contig=<ID=1>
##contig=<ID=2>
##contig=<ID=3>
##contig=<ID=4>
##contig=<ID=5>
##contig=<ID=6>
##contig=<ID=7>
##contig=<ID=8>
##contig=<ID=9>
##contig=<ID=10>
##contig=<ID=11>
##contig=<ID=12>
##contig=<ID=13>
##contig=<ID=14>
##contig=<ID=15>
##contig=<ID=16>
##contig=<ID=17>
##contig=<ID=18>
##contig=<ID=19>
##contig=<ID=20>
##contig=<ID=21>
##contig=<ID=22>
##contig=<ID=X>
##contig=<ID=Y>
##FILTER=<ID=NM,Description="Not mutation">
##FILTER=<ID=AN,Description="Either tumor or normal lacks alleles">
##FILTER=<ID=RR,Description="Reversion to reference">
##FILTER=<ID=FFT,Description="Failed Fisher Test">
##FILTER=<ID=TMA,Description="Too many alleles">
##FILTER=<ID=INS,Description="Insufficient support">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"""
    with open(output_vcf_fp, 'w') as ovf:
        ovf.write(header + "\n" + vcf_lines)


if __name__ == '__main__':
    convert_tsv_to_vcf("/home/avraham/MaruvkaLab/msmutect_runs/results/mod.partial.mut.tsv",
                       "/home/avraham/MaruvkaLab/msmutect_runs/results/mod.vcf")
