# cython: language_level=3
import scipy.stats as stats
import numpy as np

from src.IndelCalling.AlleleSet import AlleleSet
from src.IndelCalling.Histogram import Histogram


class AllelesMaximumLikelihood:
    def __init__(self, histogram: Histogram, proper_lengths: np.array, supported_lengths: np.array, noise_table: np.matrix, num_alleles: int):
        # supported lengths is a subset of proper lengths (supported must also have 6+ supporting reads)
        # gets  all lengths with at least 5 read support
        self.histogram = histogram
        self.repeat_lengths = proper_lengths
        self.num_reads = np.array([histogram.rounded_repeat_lengths[length] for length in self.repeat_lengths], dtype=np.int32)
        self.supported_repeat_lengths = supported_lengths
        self.num_alleles = num_alleles #: NEW
        # self.num_alleles = supported_lengths.size
        self.noise_table = noise_table
        self.max_log_likelihood = -1e9
        self.best_alleles = np.array([])
        self.best_frequencies = np.array([])

    def reset_intermediates(self):
        randomized_order = np.random.permutation(self.supported_repeat_lengths.size)[0:self.num_alleles]
        self.random_repeat_lengths = self.supported_repeat_lengths[randomized_order]
        self.new_frequencies = np.zeros(self.num_alleles)
        self.frequencies = np.ones(self.num_alleles) / self.num_alleles
        self.Z_i_j = np.zeros([44, self.num_alleles])
        self.new_theta = np.zeros(self.num_alleles)
        self.change = 1e6
        self.prev_log_likelihood = 1e6

    def update_ZIJ(self):
        for length in self.repeat_lengths:
            for j in range(self.num_alleles):
                self.Z_i_j[length, j] = self.noise_table[self.random_repeat_lengths[j], length] * self.frequencies[j] / np.sum(self.noise_table[self.random_repeat_lengths[:], length] * self.frequencies[:] + 1e-10)

    def estimate_new_frequencies(self):
        for k in range(self.num_alleles):
            self.new_frequencies[k] = np.sum(self.Z_i_j[self.repeat_lengths, k] * self.num_reads) / np.sum(self.num_reads)

    def maximize_new_thetas(self):
        Theta_new_temp = np.zeros(self.supported_repeat_lengths.size)
        for j in range(self.num_alleles):
            for k in range(self.supported_repeat_lengths.size):
                Test_theta = self.supported_repeat_lengths[k]
                Theta_new_temp[k] = sum(self.Z_i_j[self.repeat_lengths, j] * np.log(
                    self.noise_table[Test_theta, self.repeat_lengths] + 1e-10) * self.num_reads)
            self.new_theta[j] = self.supported_repeat_lengths[Theta_new_temp.argmax()]

    def update_frequencies(self):
        for k in range(self.num_alleles):
            self.random_repeat_lengths[k] = self.new_theta[k]
            self.frequencies[k] = self.new_frequencies[k]

    def get_log_likelihood(self) -> float:
        log_likelihood = 0
        for k in range(self.repeat_lengths.size):
            log_likelihood += self.num_reads[k] * np.log(np.sum(self.frequencies * self.noise_table[self.random_repeat_lengths, self.repeat_lengths[k]]) + 1e-10)
        return log_likelihood

    def update_guess(self) -> float:
        log_likelihood = self.get_log_likelihood()
        if log_likelihood > self.max_log_likelihood:
            self.max_log_likelihood = log_likelihood
            self.best_alleles = self.random_repeat_lengths
            self.best_frequencies = self.frequencies
        change = np.abs(self.prev_log_likelihood - log_likelihood)
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


def find_alleles(histogram: Histogram, proper_lengths: np.array, supported_repeat_lengths: np.array, noise_table: np.array, min_read_support: int = -1) -> AlleleSet:
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

    proper_motif_sizes = np.array([repeat_size for repeat_size in histogram.rounded_repeat_lengths if
                                         passes_filter(len(histogram.locus.pattern), repeat_size)])
    supported_proper_motifs = np.array([length for length in proper_motif_sizes
                                        if histogram.rounded_repeat_lengths[length]>=required_read_support])
    if supported_proper_motifs.size == 0:
        return AlleleSet(histogram, log_likelihood=-1, repeat_lengths=np.array([]), frequencies=np.array([-1]), min_read_support=required_read_support)
    elif supported_proper_motifs.size == 1:
        return AlleleSet(histogram=histogram,  log_likelihood=0, repeat_lengths=np.array(list(supported_proper_motifs)), frequencies=np.array([1]), min_read_support=required_read_support)
    else:
        return find_alleles(histogram, proper_motif_sizes, supported_proper_motifs, noise_table, required_read_support)
