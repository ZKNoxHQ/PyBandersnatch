# -*- coding: utf-8 -*-
import unittest
from field import Field
from curve import Curve


class TestCurve(unittest.TestCase):

    def set_up_curve(self):
        """Creates Bandersnatch ellptic curve.

        Test vectors generated using the file `curve_test_vectors.sage`.

        """
        F = Field(
            0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001)
        a = F(0x4247698f4e32ad45a293959b4ca17afa4a2d2317e4c6ce5023e1fd63d1b5de98)
        b = F(5)
        r = 0x1cfb69d4ca675f520cce760202687600ff8f87007419047174fd06b52876e7e1
        h = 4
        E = Curve(a, b, r, h)
        test_vectors = {}
        test_vectors['p'] = E(
            F(0x3d65d0dd26c89b6c1e111a4b7d4875a8fc60cc1a2191beba4f6c11f2fb1c2cbf), 1)
        test_vectors['q'] = E(
            F(0x4f9ace66d99bf333dd1eb464c902fcc9217372d8902583ff2e65fb3df88fdbb8), 1)
        test_vectors['p_double'] = E(
            F(0x5f6aff9f5f3a3dedf24000ef9899766ccdadc2183675dac057d3d754faf6174), 1)
        test_vectors['p_plus_q'] = E(
            F(0x6770424666db6f35101a26a4ea455b8886346935cb00298962781fa1f5f654af), 1)
        test_vectors['p_minus_q'] = E(
            F(0x45d99562094ef0bfe61ffb151ce57cf76c7c3b247179cb3ebcee99c3e0a59889), 1)
        test_vectors['φ_p'] = E(
            F(0x48850dbbc447be3c1f566275495c29ad981ed87de6818a73e1918a130a9b524d), 1)
        test_vectors['φ_minus_one_p'] = E(
            F(0x1521dba053303cfa98d4e0f358cea947c69d96f5bf8a8d3a075694b4c9b6558d), 1)
        test_vectors['k_times_p'] = E(
            F(0x35d13b10631a46276be4c289a63570de6c552376605b2881e445031c90a96a7e), 1)
        test_vectors['k1_times_p'] = E(
            F(0x367b03140f02455d6d00a7924a65635f717f6354c23ccc2d5ed2b11e4b736c9d), 1)
        test_vectors['k2_times_p'] = E(
            F(0xd6643d2acb959c8f4306b8398ddf8102dbd9fcc52c17817e1bf7943cb8a40c7), 1)
        test_vectors['k1_times_p_plus_k2_times_q'] = E(
            F(0x182a5d0dce9c87ecd9050d911058c4ec30cbc680f94e5c347da06eb068c2dd15), 1)
        test_vectors['k'] = 0x1a862619b8224e61eb24bb583c84ce04913064d37308623924c7a64fcdc9f191
        test_vectors['k1'] = -0x2286ed83a0b1545d1b7788921e40bb14
        test_vectors['k2'] = 0x6886451b4aa55294c626bb34d42e242
        test_vectors['λ'] = 0x13b4f3dc4a39a493edf849562b38c72bcfc49db970a5056ed13d21408783df05
        return E, test_vectors

    def test_j_invariant(self):
        """The j-invariant of Bandersnatch is 8000."""
        E, test_vectors = self.set_up_curve()
        self.assertEqual(E.j_inv(), 8000)

    def test_cofactor(self):
        """Checks that the multiplication by the cofactor `h` returns a prime order `r` point."""
        E, _ = self.set_up_curve()
        for i in range(10):
            p = E.random().naive_mul(E.h)
            self.assertTrue(p.is_prime_order(E.r))

    def test_is_prime_order(self):
        """Check that `p` is of prime order `r`."""
        E, test_vectors = self.set_up_curve()
        self.assertTrue(test_vectors['p'].is_prime_order(E.r))

    def test_φ_norm(self):
        """Check that φ²(P) = [-2]P."""
        E, test_vectors = self.set_up_curve()
        φ2_p = test_vectors['p'].φ().φ()
        self.assertEqual(φ2_p, test_vectors['p_double'])

    def test_φ_eigenvalue(self):
        """Check that φ acts as [λ] on the corresponding eigen-space."""
        E, test_vectors = self.set_up_curve()
        self.assertEqual(test_vectors['p'].φ(),
                         test_vectors['λ']*test_vectors['p'])

    def test_φ(self):
        """φ using test vectors."""
        E, test_vectors = self.set_up_curve()
        φ_p = test_vectors['p'].φ()
        self.assertEqual(φ_p, test_vectors['φ_p'])

    def test_φ_minus_one(self):
        """φ-1 using test vectors."""
        E, test_vectors = self.set_up_curve()
        φ_minus_one_p = test_vectors['p'].φ_minus_one()
        self.assertEqual(φ_minus_one_p, test_vectors['φ_minus_one_p'])

    def test_add(self):
        """Addition using test vectors."""
        E, test_vectors = self.set_up_curve()
        q = test_vectors['q']
        p_minus_q = test_vectors['p_minus_q']
        p_plus_q = test_vectors['p'].add(q, p_minus_q)
        self.assertEqual(p_plus_q, test_vectors['p_plus_q'])

    def test_dbl(self):
        """Doubling using test vectors."""
        E, test_vectors = self.set_up_curve()
        p_double = test_vectors['p'].dbl()
        self.assertEqual(p_double, test_vectors['p_double'])

    def test_scalar_mul(self):
        """Scalar multiplication using test vectors."""
        E, test_vectors = self.set_up_curve()
        k = test_vectors['k']
        k_times_p = test_vectors['p'].naive_mul(k)
        self.assertEqual(k_times_p, test_vectors['k_times_p'])

        k1 = test_vectors['k1']
        k1_times_p = test_vectors['p'].naive_mul(k1)
        self.assertEqual(k1_times_p, test_vectors['k1_times_p'])

        k2 = test_vectors['k2']
        k2_times_p = test_vectors['p'].naive_mul(k2)
        self.assertEqual(k2_times_p, test_vectors['k2_times_p'])

    def test_multi_scalar_mul(self):
        """Multi scalar multiplication using test vectors."""
        E, test_vectors = self.set_up_curve()
        F = E.field
        q = test_vectors['q']
        p_minus_q = test_vectors['p_minus_q']
        k1 = test_vectors['k1']
        k2 = test_vectors['k2']
        k1_p_plus_k2_q = test_vectors['p'].multi_scalar_mul(
            k1, q, k2, p_minus_q)
        self.assertEqual(
            k1_p_plus_k2_q, test_vectors['k1_times_p_plus_k2_times_q'])

    def test_glv(self):
        """GLV using test vectors."""
        E, test_vectors = self.set_up_curve()
        k = test_vectors['k']
        k_times_p = test_vectors['p'].glv(k)
        self.assertEqual(k_times_p, test_vectors['k_times_p'])

    # def test_bench(self):
    #     from time import time
    #     E, test_vectors = self.set_up_curve()
    #     # k = test_vectors['k'] # this is curiously improving even better GLV...
    #     k = 1199715452959664211345788999754554038123456678896406483911385998427296165917073  # random

    #     t = time()
    #     for i in range(50):
    #         k_times_p_glv = k*test_vectors['p']  # glv
    #     time_glv = time() - t

    #     t = time()
    #     for i in range(50):
    #         k_times_p_naive = test_vectors['p'].naive_mul(k)
    #     time_naive = time()-t
    #     print("GLV is {:.0f}% faster than a scalar multiplication.".format(
    #         time_glv/time_naive*100))
    #     self.assertEqual(k_times_p_glv, k_times_p_naive)


if __name__ == '__main__':
    unittest.main()
