# -*- coding: utf-8 -*-
import unittest
from field import Field
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
file_path = os.path.join(script_dir, 'field_test_vectors.py')


class TestField(unittest.TestCase):

    def set_up_field(self):
        """Create Bandersnatch base field.

        Test vectors generated using the file `field_test_vectors.sage`.

        """
        try:
            with open(file_path, "r") as file:
                exec(file.read(), globals())
        except FileNotFoundError as e:
            raise unittest.SkipTest(
                "The file 'field_test_vectors.py' was not found. Please generate it using `sage field_test_vectors.sage > field_test_vectors.py`.")
        return F, test_vectors  # type: ignore

    def test_random(self):
        """Two random elements are different.

        It happens with probability 1/`p`, small for large `p`."""
        F, test_vectors = self.set_up_field()
        a = F.random()
        b = F.random()
        self.assertFalse(a == b)

    def test_mul(self):
        """a*b from test vectors"""
        F, test_vectors = self.set_up_field()
        a = test_vectors['a']
        b = test_vectors['b']
        a_mul_b = a*b
        self.assertEqual(a_mul_b, test_vectors['a_mul_b'])

    def test_div(self):
        """a/b from test vectors"""
        F, test_vectors = self.set_up_field()
        a = test_vectors['a']
        b = test_vectors['b']
        a_div_b = a/b
        self.assertEqual(a_div_b, test_vectors['a_div_b'])

    def test_sqrt(self):
        """Square root test using a square"""
        F, test_vectors = self.set_up_field()
        sq = test_vectors['b']
        root = sq.sqrt()
        self.assertEqual(root*root, sq)
        self.assertTrue(
            root == test_vectors['sqrt_b'] or root == -test_vectors['sqrt_b'])

    def test_is_square(self):
        """Legendre symbol test on small squares and a non-square"""
        F, test_vectors = self.set_up_field()
        self.assertFalse(test_vectors['non_square'].is_square())
        for i in range(F.non_square.value):
            self.assertTrue(F(i).is_square())

    def run_all_test(self):
        print("Tests for Field")
        print("------")
        for method_name in dir(self):
            if method_name.startswith("test_"):
                method = getattr(self, method_name)
                if callable(method):
                    print(f"Running: {method_name}")
                    method()
        print("------")
        print()
