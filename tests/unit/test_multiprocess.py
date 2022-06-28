"""Test the code used for parallel processing."""

import unittest

from phammseqs.multiprocess import parallelize


def _square(x):
    """Function that calculates and returns the square of a number.

    :param x: a number to raise to the second power
    :type x: int or float
    """
    return x * x


def _no_retval(x):
    """Function that does nothing and returns nothing.

    :param x: a number to raise to the second power
    :type x: int or float
    """


def _second_and_fourth_powers(x):
    """Function that returns the square of a number, and its square

    :param x: a number to raise to the second and fourth powers
    :type x: int or float
    """
    return x * x, x * x * x * x


n_jobs = 50
inputs = [x for x in range(n_jobs)]
squares = [_square(x) for x in range(n_jobs)]
complex_squares = [_second_and_fourth_powers(x) for x in range(n_jobs)]


class TestParallelizeYesRetval(unittest.TestCase):
    def test_parallelize_quiet_1_cpu(self):
        """Test that parallelize returns the correct values when using
        1 CPU core, and prints nothing if verbose=False.
        """
        jobs = [(x,) for x in inputs]
        results = parallelize(jobs, 1, _square, verbose=False)
        self.assertEqual(squares, sorted(results))

    def test_parallelize_quiet_2_cpu(self):
        """Test that parallelize returns the correct values when using
        2 CPU cores, and prints nothing if verbose=False.
        """
        jobs = [(x,) for x in inputs]
        results = parallelize(jobs, 2, _square, verbose=False)
        self.assertEqual(squares, sorted(results))

    def test_parallelize_quiet_4_cpu(self):
        """Test that parallelize returns the correct values when using
        4 CPU cores, and prints nothing if verbose=False.
        """
        jobs = [(x,) for x in inputs]
        results = parallelize(jobs, 4, _square, verbose=False)
        self.assertEqual(squares, sorted(results))


class TestParallelizeNoRetval(unittest.TestCase):
    def test_parallelize_quiet_1_cpu(self):
        """Test that parallelize returns an empty list when using
        1 CPU core, and prints nothing if verbose=False.
        """
        jobs = [(x,) for x in inputs]
        results = parallelize(jobs, 1, _no_retval, verbose=False)
        self.assertEqual(list(), sorted(results))

    def test_parallelize_quiet_2_cpu(self):
        """Test that parallelize returns an empty list when using
        2 CPU cores, and prints nothing if verbose=False.
        """
        jobs = [(x,) for x in inputs]
        results = parallelize(jobs, 2, _no_retval, verbose=False)
        self.assertEqual(list(), sorted(results))

    def test_parallelize_quiet_4_cpu(self):
        """Test that parallelize returns an empty list when using
        4 CPU cores, and prints nothing if verbose=False.
        """
        jobs = [(x,) for x in inputs]
        results = parallelize(jobs, 2, _no_retval, verbose=False)
        self.assertEqual(list(), sorted(results))


class TestParallelizeComplexRetval(unittest.TestCase):
    def test_parallelize_quiet_1_cpu(self):
        """Test that parallelize returns an empty list when using
        1 CPU core, and prints nothing if verbose=False.
        """
        jobs = [(x,) for x in inputs]
        results = parallelize(jobs, 1, _second_and_fourth_powers, verbose=False)
        self.assertEqual(complex_squares, sorted(results))

    def test_parallelize_quiet_2_cpu(self):
        """Test that parallelize returns an empty list when using
        2 CPU cores, and prints nothing if verbose=False.
        """
        jobs = [(x,) for x in inputs]
        results = parallelize(jobs, 2, _second_and_fourth_powers, verbose=False)
        self.assertEqual(complex_squares, sorted(results))

    def test_parallelize_quiet_4_cpu(self):
        """Test that parallelize returns an empty list when using
        4 CPU cores, and prints nothing if verbose=False.
        """
        jobs = [(x,) for x in inputs]
        results = parallelize(jobs, 2, _second_and_fourth_powers, verbose=False)
        self.assertEqual(complex_squares, sorted(results))


if __name__ == '__main__':
    unittest.main()
