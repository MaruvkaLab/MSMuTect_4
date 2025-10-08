# cython: language_level=3
import math, sys
from collections import defaultdict
import numpy as np


class Fisher:
    def __init__(self):
        self.already_computed = defaultdict(lambda: -1)
        self.already_computed[0] = 1
        self.already_computed[1] = 1

    def factorial(self, n: int) -> int:
        if self.already_computed[n] != -1:
            return self.already_computed[n]
        elif n > sys.getrecursionlimit() - 1:  # to avoid recursion overload
            ans = math.factorial(n)
            self.already_computed[n] = ans
            return ans
        else:
            answer = n * self.factorial(n-1)
            self.already_computed[n] = answer
            return answer

    def choose(self, n: int, k: int) -> int:
        numerator = self.factorial(n)
        denominator = self.factorial(k)*self.factorial(n-k)
        return numerator // denominator

    def get_mantissa(self, n: int, num_digits: int) -> int:
        # get first prefix_length digits of n
        return int(str(n)[:num_digits])

    def big_divide(self, numerator: int, denominator: int) -> float:
        # does division for massive numbers without causing overflow error
        numerator_log = math.log(numerator)
        denominator_log = math.log(denominator)
        return math.exp(numerator_log-denominator_log)
        # numerator_power = int(math.log10(numerator))
        # denominator_power = int(math.log10(denominator))
        # numerator_mantissa_power = min(numerator_power + 1, 10)
        # numerator_mantissa = self.get_mantissa(numerator, numerator_mantissa_power)
        # denominator_mantissa_power = min(denominator_power + 1, 10)
        # denominator_mantissa = self.get_mantissa(denominator, denominator_mantissa_power)
        # quotient_mantissa = numerator_mantissa / denominator_mantissa
        # quotient = quotient_mantissa * (10 ** ((numerator_power - numerator_mantissa_power) - (denominator_power - denominator_mantissa_power)))

    def test(self, first_set: np.array, second_set: np.array):
        p_value = 1
        for i in range(first_set.size):
            # casted to int, so if number is too large for numpy int 64 bits
            p_value *= self.choose(int(first_set[i] + second_set[i]), int(first_set[i]))
        p_value = self.big_divide(p_value, self.choose(int(np.sum(first_set)+np.sum(second_set)), int(np.sum(first_set))))
        return p_value


def big_divide(numerator: int, denominator: int) -> float:
    # does division for massive numbers without causing overflow error
    numerator_log = math.log(numerator)
    denominator_log = math.log(denominator)
    return math.exp(numerator_log-denominator_log)


def one_sided_fisher_test(first_set: np.array, second_set: np.array):
    p_value = 1
    for i in range(first_set.size):
        p_value *= math.comb(int(first_set[i] + second_set[i]), int(first_set[i]))
    p_value = big_divide(p_value, math.comb(int(np.sum(first_set) + np.sum(second_set)), int(np.sum(first_set))))
    return p_value


if __name__ == '__main__':
    f=Fisher()
    a=f.big_divide(10, 3)
    print(a)