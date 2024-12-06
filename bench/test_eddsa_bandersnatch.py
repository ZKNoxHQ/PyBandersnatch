# -*- coding: utf-8 -*-
from timeit import timeit
import unittest
from src.primitives.eddsa import EdDSA


class BenchEdDSABandersnatch(unittest.TestCase):

    def set_up_signature(self):
        """Creates Bandersnatch elliptic curve.

        Test vectors generated using `sage sage/edwards.sage > tests/vectors/edwards.py`.

        """
        with open("tests/vectors/edwards-bandersnatch.py", "r") as file:
            exec(file.read(), globals())
        return EdDSA(E)  # type: ignore

    def test_bench_sign(self):
        """Benchmark signature computation."""
        global user
        user = self.set_up_signature()
        n_iter = 50
        sign_time = timeit(
            "s = user.sign(\"This is a benchmark of a signature computation\")", globals=globals(), number=n_iter)
        print("Bandersnatch - DSA\n\tSignature computation time: {:.2f}ms".format(
            sign_time/n_iter * 10**3))

    def test_bench_verif(self):
        """Benchmark signature verification."""
        global user, sig
        user = self.set_up_signature()
        sig = user.sign("This is a benchmark of a signature verification")
        n_iter = 50
        verif_time = timeit(
            "v = user.verify(\"This is a benchmark of a signature verification\", sig)", globals=globals(), number=n_iter)
        print("Bandersnatch - DSA\n\tSignature verification time: {:.2f}ms".format(
            verif_time/n_iter * 10**3))
