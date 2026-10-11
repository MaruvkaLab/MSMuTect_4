# cython: language_level=3
import math, time, random
import scipy.stats as stats
import numpy as np
from typing import List

from src.GenomicUtils.NoiseTable import get_noise_table
from src.IndelCalling.AlleleSet import AlleleSet
from src.IndelCalling.Histogram import Histogram

# import warnings
# warnings.filterwarnings("error")

random.seed(402)


class AllelesMaximumLikelihood:
    def __init__(self, histogram: Histogram, proper_lengths: List[int], supported_lengths: List[int],
                 noise_table: List[List[float]], num_alleles: int):
        # supported lengths is a subset of proper lengths (supported must also have 6+ supporting reads)
        # gets  all lengths with at least 5 read support
        self.histogram = histogram
        self.repeat_lengths = proper_lengths
        self.num_reads = [int(histogram.rounded_repeat_lengths[length]) for length in self.repeat_lengths]
        self.supported_repeat_lengths = supported_lengths
        self.num_alleles = num_alleles #: NEW
        # self.num_alleles = supported_lengths.size
        self.noise_table = noise_table
        self.max_log_likelihood = -1e9
        self.best_alleles = []
        self.best_frequencies = []

    def random_order(self, n: int, num_el: int):
        a=list(range(n))
        random.shuffle(a)
        return a[:num_el]

    def reset_intermediates(self):
        randomized_order = self.random_order(len(self.supported_repeat_lengths), self.num_alleles) #np.random.permutation(len(self.supported_repeat_lengths))[0:self.num_alleles]
        self.random_repeat_lengths = [self.supported_repeat_lengths[r] for r in randomized_order]#self.supported_repeat_lengths[randomized_order]
        self.new_frequencies = [0] * self.num_alleles
        self.frequencies = [1/self.num_alleles] * self.num_alleles
        base_list_for_Z_i_j = [0] * self.num_alleles
        self.Z_i_j: List[List[float]] = [base_list_for_Z_i_j.copy() for i in range(44)]##np.zeros([44, self.num_alleles])
        self.new_theta = [0] * self.num_alleles
        self.change = 1e6
        self.prev_log_likelihood = 1e6

    def update_ZIJ(self):
        for length in self.repeat_lengths:
            for j in range(self.num_alleles):
                new_summed = 1e-10
                i_length = int(length)
                for k in range(len(self.frequencies)):
                    new_summed+=self.frequencies[k]*self.noise_table[int(self.random_repeat_lengths[k])][i_length]
                # summed = np.sum(self.noise_table[self.random_repeat_lengths[:], length] * self.frequencies[:] + 1e-10)
                self.Z_i_j[i_length][j] = (self.noise_table[int(self.random_repeat_lengths[j])][i_length] * self.frequencies[j]) / new_summed
                # self.Z_i_j[length, j] = self.noise_table[self.random_repeat_lengths[j], length] * self.frequencies[j] / np.sum(self.noise_table[self.random_repeat_lengths[:], length] * self.frequencies[:] + 1e-10)

    def estimate_new_frequencies(self):
        for k in range(self.num_alleles):
            new_summed = 0
            for i, r in enumerate(self.repeat_lengths):
                new_summed+=self.Z_i_j[int(r)][k] * self.num_reads[i]
            self.new_frequencies[k] = new_summed/sum(self.num_reads)

            # self.new_frequencies[k] = np.sum(self.Z_i_j[self.repeat_lengths, k] * self.num_reads) / np.sum(self.num_reads)

    def maximize_new_thetas(self):
        Theta_new_temp = [0] * len(self.supported_repeat_lengths)#np.zeros(self.supported_repeat_lengths.size)
        for j in range(self.num_alleles):
            for k in range(len(self.supported_repeat_lengths)):#self.supported_repeat_lengths.size):
                Test_theta = int(self.supported_repeat_lengths[k])
                summed = 1e-10
                for i, r in enumerate(self.repeat_lengths):
                    summed+=self.Z_i_j[int(r)][j] * math.log(self.noise_table[Test_theta][int(r)]+1e-10) * self.num_reads[i]
                # Theta_new_temp[k] = sum(self.Z_i_j[self.repeat_lengths, j] * np.log(
                #     self.noise_table[Test_theta, self.repeat_lengths] + 1e-10) * self.num_reads)
                Theta_new_temp[k] = summed
            Theta_new_temp_argmax = max(range(len(Theta_new_temp)), key=Theta_new_temp.__getitem__) # thanks to SO user gg349
            self.new_theta[j] = self.supported_repeat_lengths[Theta_new_temp_argmax]

    def update_frequencies(self):
        for k in range(self.num_alleles):
            self.random_repeat_lengths[k] = self.new_theta[k]
            self.frequencies[k] = self.new_frequencies[k]

    def get_log_likelihood(self) -> float:
        log_likelihood = 0
        for k in range(len(self.repeat_lengths)): #self.repeat_lengths.size):
            summed = 1e-10
            for i,r in enumerate(self.random_repeat_lengths):
                summed+=self.frequencies[i]*self.noise_table[int(r)][int(self.repeat_lengths[k])]
            summed_log = math.log(summed)
            log_likelihood += summed_log*self.num_reads[k]
        return log_likelihood

    def update_guess(self) -> float:
        log_likelihood = self.get_log_likelihood()
        if log_likelihood > self.max_log_likelihood:
            self.max_log_likelihood = log_likelihood
            self.best_alleles = self.random_repeat_lengths
            self.best_frequencies = self.frequencies
        change = abs(self.prev_log_likelihood - log_likelihood)
        self.prev_log_likelihood = log_likelihood
        return change

    def get_alleles(self):
        for _ in range(10):
            self.reset_intermediates()
            while self.change > 1e-5: # if we have already converged, return
                self.update_ZIJ()
                self.estimate_new_frequencies()
                self.maximize_new_thetas()
                self.update_frequencies()
                self.change = self.update_guess()

        return AlleleSet(histogram=self.histogram, log_likelihood=self.max_log_likelihood,
                         repeat_lengths = self.best_alleles, frequencies=self.best_frequencies)


def find_alleles(histogram: Histogram, proper_lengths: List[int], supported_repeat_lengths: List[int], noise_table: List[List[float]], min_read_support: int = -1) -> AlleleSet:
    lesser_alleles_set = AllelesMaximumLikelihood(histogram, proper_lengths, supported_repeat_lengths, noise_table, num_alleles=1).get_alleles()
    lesser_alleles_set.min_read_support = min_read_support
    for i in range(2, 5):
        greater_alleles_set = AllelesMaximumLikelihood(histogram, proper_lengths, supported_repeat_lengths, noise_table, num_alleles=i).get_alleles()
        greater_alleles_set.min_read_support = min_read_support
        likelihood_increase = 2 * (greater_alleles_set.log_likelihood - lesser_alleles_set.log_likelihood)
        if likelihood_increase > 0:
            p_value_i_alleles = stats.chi2.pdf(likelihood_increase, 2)
            if p_value_i_alleles > 0.05:
                return lesser_alleles_set
            elif len(supported_repeat_lengths) == i:
                return greater_alleles_set
            else:
                lesser_alleles_set = greater_alleles_set
        else:  # should only ever occur on first time through (ie. 1 and 2 allele sets)
            return lesser_alleles_set
    return lesser_alleles_set


def repeat_threshold(ms_length: int):
    # number of repeats necessary for microsatellite of given length to be considered
    if ms_length == 1:
        return 5
    elif ms_length == 2:
        return 4
    elif ms_length >= 3:
        return 3


def passes_filter(motif_length: int, repeat_size: float):
    return repeat_threshold(motif_length) <= repeat_size <= 40


def calculate_alleles(histogram: Histogram, noise_table: np.array, required_read_support):
    if len(histogram.repeat_lengths)==0:
        return AlleleSet(histogram, log_likelihood=-1, repeat_lengths=np.array([]), frequencies=np.array([-1]), min_read_support=required_read_support)

    proper_motif_sizes = [repeat_size for repeat_size in histogram.rounded_repeat_lengths if
                                         passes_filter(len(histogram.locus.pattern), repeat_size)]
    supported_proper_motifs = [length for length in proper_motif_sizes
                                        if histogram.rounded_repeat_lengths[length]>=required_read_support]
    if len(supported_proper_motifs) == 0:
        return AlleleSet(histogram, log_likelihood=-1, repeat_lengths=np.array([]), frequencies=np.array([-1]), min_read_support=required_read_support)
    elif len(supported_proper_motifs) == 1:
        return AlleleSet(histogram=histogram,  log_likelihood=0, repeat_lengths=np.array(list(supported_proper_motifs)), frequencies=np.array([1]), min_read_support=required_read_support)
    else:
        ret = find_alleles(histogram, proper_motif_sizes, supported_proper_motifs, noise_table, required_read_support)
        ret.frequencies = np.array(ret.frequencies) # to keep type consistency
        ret.repeat_lengths = np.array(ret.repeat_lengths)
        return ret


if __name__ == '__main__':
    nt = get_noise_table()
    td_nt = nt.tolist()
    random_indices = [(random.randint(0, 40), random.randint(0, 40)) for i in range(1_000_000)]
    np_test_start = time.time()
    for i in range(1_000_000):
        randindx = random_indices[i]
        _ = nt[randindx[0], randindx[1]]
    np_test_end = time.time()
    print(f"NP TEST took {np_test_end-np_test_start}")
    print(type(nt))
    np_test_start = time.time()
    for i in range(1_000_000):
        randindx = random_indices[i]
        _ = td_nt[randindx[0]][randindx[1]]
    np_test_end = time.time()
    print(f"NP TEST took {np_test_end - np_test_start}")