# -*- coding: utf-8 -*-
from timeit import timeit
import unittest


class BenchEdwards25519(unittest.TestCase):

    def set_up_curve(self):
        """Creates Bandersnatch elliptic curve.

        Test vectors generated using `sage sage/edwards.sage > tests/vectors/edwards.py`.

        """
        with open("tests/vectors/edwards-25519.py", "r") as file:
            exec(file.read(), globals())
        return E, test_vectors  # type: ignore

    def test_bench_multi_scalar_mul_2(self):
        """Benchmark GLV and the naive scalar multiplication using `p` and `k`."""
        E, test_vectors = self.set_up_curve()
        global k, p, q
        k = test_vectors['k']
        p = test_vectors['p']
        q = test_vectors['q']

        n_iter = 50
        naive_mul_time = timeit("res_1 = (k*p).add(k*q)",
                                globals=globals(), number=n_iter)

        multi_scalar_mul_time = timeit("res_2 = p.multi_scalar_mul_2(k,q,k)",
                                       globals=globals(), number=n_iter)

        print("Edwards curve 25519:\n\tNaive multi scalar mul: {:.2f} ms;\n\tMSM: {:.2f} ms ({:.0f}% faster)".format(
            naive_mul_time/n_iter*10**3, multi_scalar_mul_time / n_iter*10**3, (naive_mul_time - multi_scalar_mul_time)/naive_mul_time*100))

    def test_bench_glv(self):
        """Benchmark GLV and the naive scalar multiplication using `p` and `k`."""
        E, test_vectors = self.set_up_curve()
        k = test_vectors['k']

        n_iter = 50
        naive_mul_time = timeit("test_vectors['p'].naive_mul(test_vectors['k'])",
                                globals=globals(), number=n_iter)

        glv_time = timeit("test_vectors['k']*test_vectors['p']",
                          globals=globals(), number=n_iter)

        print("Edwards curve 25519:\n\tNaive mul: {:.2f} ms;\n\tGLV: {:.2f} ms ({:.0f}% faster)".format(
            naive_mul_time/n_iter*10**3, glv_time / n_iter*10**3, (naive_mul_time-glv_time)/naive_mul_time*100))
