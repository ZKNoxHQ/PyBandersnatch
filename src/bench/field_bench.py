# -*- coding: utf-8 -*-
import unittest
from field import Field
import sys
from random import randint
from timeit import timeit


class BenchField:

    def set_up_field(self):
        """Create Bandersnatch base field.

        Test vectors generated using the file `field_test_vectors.sage`.

        """
        try:
            with open("src/tests/field_test_vectors.py", "r") as file:
                exec(file.read(), globals())
        except FileNotFoundError as e:
            raise unittest.SkipTest(
                "The file 'field_test_vectors.py' was not found. Please generate it using `sage field_test_vectors.sage > field_test_vectors.py`.")
        return F, test_vectors  # type: ignore

    def bench_arithmetic(self):
        """Benchmark field multiplication."""
        F, test_vectors = self.set_up_field()
        n_iter = 50
        naive_mul_time = timeit(
            "test_vectors['a']*test_vectors['b']", globals=globals(), number=n_iter)
        naive_add_time = timeit(
            "test_vectors['a']+test_vectors['b']", globals=globals(), number=n_iter)
        naive_div_time = timeit(
            "test_vectors['a']/test_vectors['b']", globals=globals(), number=n_iter)

        print("Mul: {:.3f}μs".format(naive_mul_time/n_iter * 10**6))
        print("Add: {:.3f}μs".format(naive_add_time/n_iter * 10**6))
        print("Div: {:.3f}μs".format(naive_div_time/n_iter * 10**6))

    def run_all_bench(self):
        for method_name in dir(self):
            if method_name.startswith("bench_"):
                method = getattr(self, method_name)
                if callable(method):
                    print(f"Running: {method_name}")
                    method()
