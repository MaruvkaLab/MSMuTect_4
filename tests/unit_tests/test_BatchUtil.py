import unittest

from src.Entry.BatchUtil import *


class TestBatchUtil(unittest.TestCase):

    def test_BatchSize(self):
        self.assertEqual(get_batch_sizes(100_001, 100_000)[0], 100_000)
        self.assertEqual(get_batch_sizes(100_001, 100_000)[1], 1)


if __name__ == '__main__':
    unittest.main()
