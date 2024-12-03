# -*- coding: utf-8 -*-
import unittest
from src.field import Field
from src.curve.edwards import Edwards
from random import randint


class TestEdwards(unittest.TestCase):

    def set_up_curve(self):
        """Creates Bandersnatch elliptic curve.

        Test vectors generated using the file `tests/vectors/edwards.sage`.

        """
        try:
            with open('tests/vectors/edwards.py', "r") as file:
                exec(file.read(), globals())
        except FileNotFoundError as e:
            raise unittest.SkipTest(
                "The file 'tests/vectors/edwards.py' was not found. Please generate it using `sage Edwards_edwards_test_vectors.sage > Edwards_edwards_test_vectors.py`.")
        return E, test_vectors  # type: ignore

    def test_j_invariant(self):
        """j(Bandersnatch) = 8000"""
        E, test_vectors = self.set_up_curve()
        self.assertEqual(E.j_inv(), 8000)

    def test_in_curve(self):
        """Point is on the curve"""
        E, test_vectors = self.set_up_curve()
        self.assertTrue(test_vectors['p'].in_curve())
        self.assertTrue(test_vectors['q'].in_curve())
        self.assertTrue(test_vectors['p_plus_q'].in_curve())
        self.assertTrue(test_vectors['φ_p'].in_curve())
        self.assertTrue(test_vectors['k1_times_p'].in_curve())
        self.assertTrue(test_vectors['k2_times_p'].in_curve())
        self.assertTrue(test_vectors['k_times_p'].in_curve())

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

    def test_φ_norm(self):
        """φ²(p) = [-2]p"""
        E, test_vectors = self.set_up_curve()
        φ2_p = test_vectors['p'].φ().φ()
        self.assertEqual(φ2_p, test_vectors['p'].dbl().neg())

    def test_φ_eigenvalue(self):
        """φ=[λ] on the <p> eigen-space"""
        E, test_vectors = self.set_up_curve()
        self.assertEqual(test_vectors['p'].φ(),
                         test_vectors['p'].naive_mul(test_vectors['λ']))

    def test_φ(self):
        """φ from test vectors"""
        E, test_vectors = self.set_up_curve()
        φ_p = test_vectors['p'].φ()
        self.assertEqual(φ_p, test_vectors['φ_p'])

    def test_add(self):
        """p+q from test vectors"""
        E, test_vectors = self.set_up_curve()
        q = test_vectors['q']
        p_plus_q = test_vectors['p'].add(q)
        self.assertEqual(p_plus_q, test_vectors['p_plus_q'])

    def test_dbl(self):
        """2*p from test vectors"""
        E, test_vectors = self.set_up_curve()
        p_double = test_vectors['p'].dbl()
        self.assertEqual(p_double, test_vectors['p_double'])

    def test_dbl_add(self):
        """p.dbl() == p.add(p)"""
        E, test_vectors = self.set_up_curve()
        p = test_vectors['p']
        self.assertEqual(p.dbl(), p.add(p))

    def test_neg(self):
        """p + (-p) = 0"""
        E, test_vectors = self.set_up_curve()
        p = test_vectors['p']
        minus_p = p.neg()
        self.assertEqual(p.add(minus_p), E(0, 1, 1))

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

    def test_multi_scalar_mul(self):
        """k1*p + k2*q from test vectors"""
        E, test_vectors = self.set_up_curve()
        q = test_vectors['q']
        k1 = test_vectors['k1']
        k2 = test_vectors['k2']
        k1_p_plus_k2_q = test_vectors['p'].multi_scalar_mul(
            k1, q, k2)
        self.assertEqual(
            k1_p_plus_k2_q, test_vectors['k1_times_p_plus_k2_times_q'])

    def test_multi_scalar_mul_small_scalars(self):
        """k1*p + k2*q using small scalars -3 ≤ k1,k2 < 3"""
        E, test_vectors = self.set_up_curve()
        p = test_vectors['p']
        q = test_vectors['q']
        for k1 in range(-3, 3):
            for k2 in range(-3, 3):
                tmp1 = p.multi_scalar_mul(k1, q, k2)
                self.assertEqual(tmp1, p.naive_mul(k1).add(q.naive_mul(k2)))

    def test_glv(self):
        """GLV technique for k*p from test vectors"""
        E, test_vectors = self.set_up_curve()
        k = test_vectors['k']
        k_times_p = test_vectors['p'].glv(k)
        self.assertEqual(k_times_p, test_vectors['k_times_p'])

    def test_scalar_mul_random(self):
        """k1*p + k2*q from random k1,k2"""
        E, test_vectors = self.set_up_curve()
        for i in range(30):
            k = randint(1, 1 << 254)
            k_times_p_1 = test_vectors['p'].glv(k)
            k_times_p_2 = test_vectors['p'].naive_mul(k)
            self.assertEqual(k_times_p_1, k_times_p_2)

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
        self.assertEqual(-12345*test_vectors['p'],
                         12345*test_vectors['p'].neg())

    def test_generator(self):
        """Generation of the small x order r point as in the test vectors"""
        E, test_vectors = self.set_up_curve()
        # Small x generator
        x = 0
        g = E(E.field(0), E.field(1), E.field(1))
        while not (g.is_prime_order(E.r)):
            x += 1
            while not ((1-E.a*x**2)/(1-E.d*x**2)).is_square():
                x = -x
                if x > 0:
                    x += 1
            y = ((1-E.a*x**2)/(1-E.d*x**2)).sqrt()
            # TODO: what about -y?
            g = E(E.field(x), E.field(y.value), E.field(1))
        self.assertTrue(g.is_prime_order(E.r))
        self.assertEqual(E.generator, g)

    def test_small_x_point(self):
        """Point with x=1 is not of order r, from the test vectors"""
        E, test_vectors = self.set_up_curve()
        small_p = test_vectors['small_p']
        small_p_dbl = test_vectors['small_p_dbl']
        self.assertEqual(small_p.dbl(), small_p_dbl)
        self.assertEqual(small_p.naive_mul(E.r), E(0, -1, 1))
        self.assertTrue(E.generator.is_prime_order(E.r))

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
