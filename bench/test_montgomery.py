# -*- coding: utf-8 -*-
from timeit import timeit
import unittest


class BenchMontgomery(unittest.TestCase):

    def set_up_curve(self):
        """Creates Bandersnatch elliptic curve.

        Test vectors generated using `sage sage/montgomery.sage > tests/vectors/montgomery.py`.

        """
        with open("tests/vectors/montgomery.py", "r") as file:
            exec(file.read(), globals())
        return E, test_vectors  # type: ignore

    def test_bench_glv(self):
        """Benchmark GLV and the naive scalar multiplication using `p` and `k`."""
        E, test_vectors = self.set_up_curve()
        k = test_vectors['k']

        n_iter = 50
        naive_mul_time = timeit("test_vectors['p'].naive_mul(test_vectors['k'])",
                                globals=globals(), number=n_iter)

        glv_time = timeit("test_vectors['k']*test_vectors['p']",
                          globals=globals(), number=n_iter)

        print("Naive mul: {:.2f} ms; GLV: {:.2f} ms ({:.0f}% faster)".format(
            naive_mul_time/n_iter*10**3, glv_time / n_iter*10**3, glv_time/naive_mul_time*100))
