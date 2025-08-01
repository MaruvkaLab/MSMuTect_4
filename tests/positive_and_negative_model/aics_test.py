import random
from scipy.stats import binom, norminvgauss
import matplotlib.pyplot as plt
from scipy.stats import poisson

from src.GenomicUtils.NoiseTable import get_noise_table
from src.IndelCalling.CallAllelesFast import calculate_alleles
from src.IndelCalling.CallMutations import calculate_AICs, call_verified_locus, passes_AICs, fisher_test
from src.IndelCalling.FisherTest import Fisher
from tests.testing_utils.generate_histograms import get_mutation_histograms, get_mutation_histograms_parameters
import numpy as np

from tests.testing_utils.read_results import LocusMutationCall


def main():
    noise_table = get_noise_table()
    # tumor, normal = get_mutation_histograms()
    num_tumor_reads = 50
    aiks = []
    fishies = []
    normal = get_mutation_histograms_parameters([11], [30], 1)
    normal_alleles = calculate_alleles(normal, noise_table, 5)
    maximum_mut_strength = int((3/4)*num_tumor_reads)
    full_size_aiks = np.zeros((maximum_mut_strength, maximum_mut_strength))
    full_size_fish = np.zeros((maximum_mut_strength, maximum_mut_strength))
    for i in np.arange(1, maximum_mut_strength)/num_tumor_reads:
        for j in np.arange(1/num_tumor_reads, i/2, 1/num_tumor_reads):

            tumor = get_mutation_histograms_parameters([11, 10, 9],[int((1 - i) * num_tumor_reads), max(int(num_tumor_reads * (i-j)), 1),
                                                                             max(num_tumor_reads*j, 1)], 1)
            tumor_alleles = calculate_alleles(tumor, noise_table, 5)
            aik, fish = passes_AICs_and_fisher(normal_alleles, tumor_alleles, noise_table, Fisher())
            if (i>0.3) and not aik:
                croc=2+2
                fake = calculate_alleles(tumor, noise_table, 5)
            full_size_aiks[int((i)*num_tumor_reads), int((j)*num_tumor_reads)] = aik
            full_size_fish[int((i)*num_tumor_reads), int((j)*num_tumor_reads)] = fish

    # plt.plot(np.arange(0.01, 0.6, 0.01), fishies, label="fishies")
    # plt.plot(np.arange(0.01, 0.6, 0.01), aiks, label="aiks")
    # plt.legend()
    plt.imshow(full_size_aiks)
    plt.title("AIC")
    plt.show()
    plt.imshow(full_size_fish)
    plt.title("FISH")
    plt.show()
    # # print(tumor_alleles.sorted_alleles())
    # # print(call_verified_locus(normal_alleles, tumor_alleles, noise_table, Fisher()))





def main_2():
    noise_table = get_noise_table()
    # tumor, normal = get_mutation_histograms()
    num_tumor_reads = 40
    normal = get_mutation_histograms_parameters([11], [30])
    normal_alleles = calculate_alleles(normal, noise_table, 5)
    tumor = get_mutation_histograms_parameters([11, 10, 9, 6, 5], [70, 4, 4, 4, 4])
    tumor_alleles = calculate_alleles(tumor, noise_table, 5)
    aik, fish = passes_AICs_and_fisher(normal_alleles, tumor_alleles, noise_table, Fisher())
    print(aik, fish)

def passes_AICs_and_fisher(normal_alleles, tumor_alleles, noise_table, fisher_calculator):
    return passes_AICs(calculate_AICs(normal_alleles, tumor_alleles, noise_table), 8.0), fisher_test(normal_alleles, tumor_alleles, fisher_calculator) < 0.031


def median_read_coverage(poisson_mu):

    median = poisson.ppf(0.5, poisson_mu)
    return median


def main_fixed():
    noise_table = get_noise_table()
    # tumor, normal = get_mutation_histograms()
    num_tumor_reads = 95
    aiks = []
    fishies = []
    base_num_repeats = 7
    normal = get_mutation_histograms_parameters([base_num_repeats], [num_tumor_reads//3])
    normal_alleles = calculate_alleles(normal, noise_table, 5)
    maximum_mut_strength = num_tumor_reads//2
    full_size_aiks = np.ones((maximum_mut_strength, maximum_mut_strength))*0.5
    full_size_fish = np.ones((maximum_mut_strength, maximum_mut_strength))*0.5
    for i in range(maximum_mut_strength):
        for j in range(maximum_mut_strength):

            tumor = get_mutation_histograms_parameters([base_num_repeats, base_num_repeats-1, base_num_repeats-2],
                                                       [num_tumor_reads-i-j, max(i-j, 1), max(j, 1)])
            tumor_alleles = calculate_alleles(tumor, noise_table, 5)
            aik, fish = passes_AICs_and_fisher(normal_alleles, tumor_alleles, noise_table, Fisher())
            if (j>40) and not aik:
                croc=2+2
                calculate_alleles(tumor, noise_table, 5)
            full_size_aiks[i, j] = int(aik)
            full_size_fish[i, j] = int(fish)

    plt.imshow(full_size_aiks)
    plt.title("AIC")
    plt.xlabel("Larger Mutation Reads")
    plt.ylabel("Smaller Mutation Reads")
    plt.show()
    plt.imshow(full_size_fish)
    plt.xlabel("Larger Mutation Reads")
    plt.ylabel("Smaller Mutation Reads")
    plt.title("FISHER")
    plt.show()


def neccesary_read_coverage(noise_table, base_num_repeats, new_num_repeats, num_reads, normal_alleles, pattern_length) -> int:
    mutation=False
    num_mut_reads = 5
    while num_mut_reads<num_reads//2+1:
        tumor = get_mutation_histograms_parameters([base_num_repeats, new_num_repeats],
                                               [num_reads-num_mut_reads, num_mut_reads], pattern_length)
        tumor_alleles = calculate_alleles(tumor, noise_table, 5)

        aik, fish = passes_AICs_and_fisher(normal_alleles, tumor_alleles, noise_table, Fisher())
        if len(tumor_alleles.repeat_lengths)>1:
            croc=1
        mutation = aik and fish
        if mutation:
            return num_mut_reads
        else:
            num_mut_reads+=1
    return 1

def td_graph():
    median_locus_lengths = [int(c) for c in [5.0, 11.0, 11.0, 15.0, 17.0, 21.0, 25.0, 29.0, 33.0, 36.0, 35.0, 45.0, 44.0, 50.0, 60.0]]
    depth = 80
    read_length = 151
    flanking = 10
    repeat_lengths = list(range(1, 16))
    full_size_calls = np.ones((15, 15))
    noise_table = get_noise_table()
    for repeat_length in repeat_lengths:
        median_locus_length = median_locus_lengths[repeat_length-1]
        possible_bases = (read_length - median_locus_length - 2*flanking)
        # print(possible_bases)
        poisson_mu = (depth / read_length) * possible_bases
        # print(poisson_mu)
        median_read_depth = median_read_coverage(poisson_mu)
        print(median_read_depth)
        # print("*****************")
        base_num_repeats = median_locus_length//repeat_length
        # print(base_num_repeats)
        normal = get_mutation_histograms_parameters([base_num_repeats], [median_read_depth // 3], repeat_length)
        normal_alleles = calculate_alleles(normal, noise_table, 5)
        for j in range(1, 16):
            if j == repeat_length:
                continue
            if repeat_length==6 and j == 14:
                croc=1
            req_read_coverage = neccesary_read_coverage(noise_table, base_num_repeats, j, median_read_depth, normal_alleles, repeat_length)
            full_size_calls[repeat_length-1, j-1] = req_read_coverage
            # print(repeat_length-1)
            # print(j-1)
    plt.imshow(full_size_calls)
    plt.show()

def neccesary_read_coverage_fixed(noise_table, base_num_repeats, new_num_repeats, num_reads) -> int:
    normal = get_mutation_histograms_parameters([base_num_repeats], [num_reads // 3], pattern_length=5)
    normal_alleles = calculate_alleles(normal, noise_table, 5)
    mutation=False
    num_mut_reads = 5
    aik=fish=False
    while num_mut_reads<num_reads//2+1:
        tumor = get_mutation_histograms_parameters([base_num_repeats, new_num_repeats],
                                               [num_reads-num_mut_reads, num_mut_reads], pattern_length=5)
        tumor_alleles = calculate_alleles(tumor, noise_table, 5)
        aik, fish = passes_AICs_and_fisher(normal_alleles, tumor_alleles, noise_table, Fisher())
        if base_num_repeats==6 and new_num_repeats==30:
            aik, fish = passes_AICs_and_fisher(normal_alleles, tumor_alleles, noise_table, Fisher())
        if len(tumor_alleles.repeat_lengths)>1:
            croc=1
        mutation = fish
        if mutation:
            return num_mut_reads
        else:
            num_mut_reads+=1
    return 1


def td_graph_fixed():
    median_locus_lengths = [int(c) for c in [5.0, 11.0, 11.0, 15.0, 17.0, 21.0, 25.0, 29.0, 33.0, 36.0, 35.0, 45.0, 44.0, 50.0, 60.0]]
    depth = 80
    read_length = 151
    flanking = 10
    repeat_lengths = list(range(1, 16))
    full_size_calls = np.ones((40, 40))
    noise_table = get_noise_table()
    plt.imshow(noise_table)
    plt.colorbar()
    plt.show()
    for repeat_length in repeat_lengths:
        median_locus_length = median_locus_lengths[repeat_length-1]
        possible_bases = (read_length - median_locus_length - 2*flanking)
        # print(possible_bases)
        poisson_mu = (depth / read_length) * possible_bases
        # print(poisson_mu)
        median_read_depth = median_read_coverage(poisson_mu)
        print(median_read_depth)
        # print("*****************")
        base_num_repeats = median_locus_length//repeat_length
        # print(base_num_repeats)
        for i in range(2, 40):
            for j in range(2, 40):
                if j == i:
                    continue
                # if repeat_length==6 and j == 14:
                #     croc=1
                req_read_coverage = neccesary_read_coverage_fixed(noise_table, i, j, median_read_depth)
                full_size_calls[i-1, j-1] = req_read_coverage
                # print(repeat_length-1)
                # print(j-1)
        plt.imshow(full_size_calls)
        plt.colorbar()
        plt.show()


td_graph_fixed()