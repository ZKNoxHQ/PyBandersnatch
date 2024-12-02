# -*- coding: utf-8 -*-
from timeit import timeit
import unittest
from src.primitives.eddsa import EdDSA


class BenchEdDSA(unittest.TestCase):

    def set_up_signature(self):
        """Creates Bandersnatch elliptic curve.

        Test vectors generated using `sage sage/edwards.sage > tests/vectors/edwards.py`.

        """
        with open("tests/vectors/edwards.py", "r") as file:
            exec(file.read(), globals())
        return EdDSA(E)  # type: ignore

    def test_bench_eddsa_sign(self):
        """Benchmark GLV and the naive scalar multiplication using `p` and `k`."""
        global user
        user = self.set_up_signature()
        n_iter = 50
        sign_time = timeit(
            "s = user.sign(\"This is a benchmark of a signature computation\")", globals=globals(), number=n_iter)
        print("Signature time: {:.2f}ms".format(sign_time/n_iter * 10**3))
