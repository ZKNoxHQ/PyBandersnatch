# -*- coding: utf-8 -*-
from field import Field
from curve import Curve
import sys
from random import randint
from timeit import timeit


class BenchCurve:

    def set_up_curve(self):
        """Creates Bandersnatch elliptic curve.

        Test vectors generated using the file `curve_test_vectors.sage`.

        """
        with open("curve_test_vectors.py", "r") as file:
            exec(file.read(), globals())
        return E, test_vectors  # type: ignore

    def bench_glv(self):
        """Benchmark GLV and the naive scalar multiplication using `p` and `k`."""
        E, test_vectors = self.set_up_curve()
        k = test_vectors['k']

        n_iter = 50
        naive_mul_time = timeit("test_vectors['p'].naive_mul(test_vectors['k'])",
                                globals=globals(), number=n_iter)

        glv_time = timeit("test_vectors['k']*test_vectors['p']",
                          globals=globals(), number=n_iter)

        print(f"Naive mul:{naive_mul_time/n_iter*10**6:.0f} μs ; GLV: {glv_time /
              n_iter*10**6:.0f} μs ({glv_time/naive_mul_time*100:.0f}% faster)")

    def run_all_bench(self):
        for method_name in dir(self):
            if method_name.startswith("bench_"):
                method = getattr(self, method_name)
                if callable(method):
                    print(f"Running: {method_name}")
                    method()
