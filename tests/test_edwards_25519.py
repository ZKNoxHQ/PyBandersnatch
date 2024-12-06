# -*- coding: utf-8 -*-
import unittest
from src.field import Field
from src.curve.edwards import Edwards
from random import randint


class TestEdwards25519(unittest.TestCase):

    def set_up_curve(self):
        """Creates Bandersnatch elliptic curve.

        Test vectors generated using the file `tests/vectors/edwards-25519.sage`.

        """
        try:
            with open('tests/vectors/edwards-25519.py', "r") as file:
                exec(file.read(), globals())
        except FileNotFoundError as e:
            raise unittest.SkipTest(
                "The file 'tests/vectors/edwards.py' was not found. Please generate it using `sage Edwards_edwards_test_vectors.sage > Edwards_edwards_test_vectors.py`.")
        return E, test_vectors  # type: ignore

    def test_in_curve(self):
        """Point is on the curve"""
        E, test_vectors = self.set_up_curve()
        self.assertTrue(test_vectors['p'].in_curve())
        self.assertTrue(test_vectors['q'].in_curve())
        self.assertTrue(test_vectors['p_plus_q'].in_curve())
        self.assertTrue(test_vectors['k1_times_p'].in_curve())
        self.assertTrue(test_vectors['k2_times_p'].in_curve())
        self.assertTrue(test_vectors['k_times_p'].in_curve())
        self.assertTrue(E.g.in_curve())

    def test_cofactor(self):
        """h*p is of order r"""
        E, _ = self.set_up_curve()
        for i in range(10):
            p = E.random().naive_mul(E.h)
            self.assertTrue(p.is_prime_order(E.r))

    def test_is_prime_order(self):
        """p is of prime order r"""
        E, test_vectors = self.set_up_curve()
        self.assertTrue(test_vectors['p'].is_prime_order(E.r))

    def test_add(self):
        """p+q from test vectors"""
        E, test_vectors = self.set_up_curve()
        q = test_vectors['q']
        p_plus_q = test_vectors['p'] + q
        self.assertEqual(p_plus_q, test_vectors['p_plus_q'])

    def test_dbl(self):
        """2*p from test vectors"""
        E, test_vectors = self.set_up_curve()
        p_double = test_vectors['p'].dbl()
        self.assertEqual(p_double, test_vectors['p_double'])

    def test_dbl_add(self):
        """p.dbl() == p + p"""
        E, test_vectors = self.set_up_curve()
        p = test_vectors['p']
        self.assertEqual(p.dbl(), p+p)

    def test_neg(self):
        """p + (-p) = 0"""
        E, test_vectors = self.set_up_curve()
        p = test_vectors['p']
        minus_p = p.neg()
        self.assertEqual(p + minus_p, E(0, 1, 1))

    def test_scalar_mul(self):
        """k*p from test vectors"""
        E, test_vectors = self.set_up_curve()

        k = test_vectors['k']
        p = test_vectors['p']
        k_times_p = test_vectors['p'].naive_mul(k)

        self.assertEqual(k_times_p, test_vectors['k_times_p'])

        k1 = test_vectors['k1']
        k1_times_p = test_vectors['p'].naive_mul(k1)
        self.assertEqual(k1_times_p, test_vectors['k1_times_p'])

        k2 = test_vectors['k2']
        k2_times_p = test_vectors['p'].naive_mul(k2)
        self.assertEqual(k2_times_p, test_vectors['k2_times_p'])

    def test_scalar_mul_edge_case(self):
        """α*p for {α, α-r} small"""
        E, test_vectors = self.set_up_curve()
        self.assertEqual(-2*test_vectors['p'], test_vectors['p_double'].neg())
        self.assertEqual(-1*test_vectors['p'], test_vectors['p'].neg())
        self.assertEqual(0*test_vectors['p'], E(0, 1, 1))
        self.assertEqual(1*test_vectors['p'], test_vectors['p'])
        self.assertEqual(2*test_vectors['p'], test_vectors['p_double'])
        for i in range(10):
            self.assertEqual(
                ((1 << 256)-i)*test_vectors['p'], (2**256-i-E.r) * test_vectors['p'])
            self.assertEqual((E.r+i)*test_vectors['p'], i*test_vectors['p'])

    def test_scalar_mul_negation(self):
        """k*p and -k*p"""
        E, test_vectors = self.set_up_curve()
        p = test_vectors['p']
        self.assertEqual(-12345*p,
                         12345*p.neg())

    # def test_generator(self):
    #     """Generation of the small x order r point as in the test vectors"""
    #     E, test_vectors = self.set_up_curve()
    #     # # Small x generator
    #     # x = 0
    #     # g = E(E.field(0), E.field(1), E.field(1))
    #     # while not (g.is_prime_order(E.r)):
    #     #     x += 1
    #     #     while not ((1-E.a*x**2)/(1-E.d*x**2)).is_square():
    #     #         x = -x
    #     #         if x > 0:
    #     #             x += 1
    #     #     y = ((1-E.a*x**2)/(1-E.d*x**2)).sqrt()
    #     #     # TODO: what about -y?
    #     #     g = E(E.field(x), E.field(y.value), E.field(1))
    #     # self.assertTrue(g.is_prime_order(E.r))
    #     self.assertEqual(E.generator(), g) # todo define g in test_vectors

    def test_order_2_points(self):
        """The point (0,-1,1) is of order 2"""
        E, test_vectors = self.set_up_curve()
        p_order_2_1 = test_vectors['p_order_2_1']
        self.assertTrue(p_order_2_1.is_prime_order(2))
        # TODO: We cannot manipulate points with z=0 in this implementation.
        # Hence, this does not work.
        # p_order_2_2 = test_vectors['p_order_2_2']
        # self.assertTrue(p_order_2_2.is_prime_order(2))
        # p_order_2_3 = test_vectors['p_order_2_3']
        # self.assertTrue(p_order_2_3.is_prime_order(2))

    def test_encode_decode(self):
        E, test_vectors = self.set_up_curve()
        p = test_vectors['p']
        enc_p = p.encode_base(256)
        assert E.decode_base(enc_p, 256) == p
