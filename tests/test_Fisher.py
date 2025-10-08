import unittest, math, time
import numpy as np

from src.IndelCalling.FisherTest import Fisher, one_sided_fisher_test


class TestFisher(unittest.TestCase):

    # def test_factorial(self):
    #     fisher = Fisher()
    #     self.assertEqual(int(math.log(fisher.factorial(55), 10)), 73)
    #     self.assertEqual(fisher.factorial(6), 720)
    #
    # def factorial_time(self):
    #     fisher = Fisher()
    #     start = time.process_time()
    #     for i in range(1, 100, -1):
    #         _ = fisher.factorial(i)
    #     self.assertLess(time.process_time() - start, 1e-5)
    #
    # def test_choose(self):
    #     fisher = Fisher()
    #     self.assertEqual(fisher.choose(10, 4), 210)
    #     self.assertEqual(fisher.choose(10, 7), 120)
    #
    #
    #
    # def test_size_warning(self):
    #     fisher = Fisher()
    #     a = fisher.test(np.array([100, 35, 95]), np.array([4, 6, 50]))
    #
    # def test_big_divide(self):
    #     fisher = Fisher()
    #     self.assertAlmostEqual(fisher.big_divide(200_000, 1_000), 200)
    #     self.assertAlmostEqual(fisher.big_divide(1_000, 200_000), .005)
    #     self.assertAlmostEqual(fisher.big_divide(2**3000, 2**3002), .25)


    def test_fisher_test(self):
        fisher = Fisher()
        self.assertAlmostEqual(fisher.test(np.array([1, 11], dtype=object), np.array([9, 3], dtype=object)), 0.001346076)
        self.assertAlmostEqual(one_sided_fisher_test(np.array([1, 11], dtype=object), np.array([9, 3], dtype=object)), 0.001346076)



if __name__ == '__main__':
    unittest.main()
