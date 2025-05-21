import random
from scipy.stats import binom, norminvgauss
import matplotlib.pyplot as plt

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
    normal = get_mutation_histograms_parameters([11], [30])
    normal_alleles = calculate_alleles(normal, noise_table, 5)
    maximum_mut_strength = int((3/4)*num_tumor_reads)
    full_size_aiks = np.zeros((maximum_mut_strength, maximum_mut_strength))
    full_size_fish = np.zeros((maximum_mut_strength, maximum_mut_strength))
    for i in np.arange(1, maximum_mut_strength)/num_tumor_reads:
        for j in np.arange(1/num_tumor_reads, i/2, 1/num_tumor_reads):

            tumor = get_mutation_histograms_parameters([11, 10, 9],[int((1 - i) * num_tumor_reads), max(int(num_tumor_reads * (i-j)), 1),
                                                                             max(num_tumor_reads*j, 1)])
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

main_fixed()