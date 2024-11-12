# -*- coding: utf-8 -*-
import unittest
from field import Field


class TestField(unittest.TestCase):

    def set_up_field(self):
        """Create bandersnatch base field.

        Test vectors generated using the file `field_test_vectors.sage`.

        """
        p = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
        F = Field(p)
        test_vectors = {}
        test_vectors['a'] = F(
            0x323f6ab2be9df186b32d48320a504c8841b51474108995013ca91c0b6dfaeb04)
        test_vectors['b'] = F(
            0x49f3b288d4250a92c8a6483fc259aed1cfd09151814f598f415b3f41c93c2850)
        test_vectors['a_plus_b'] = F(
            0x84575e869257ed14899b869c3082354bdc801c291da92917e045b4e37371353)
        test_vectors['a_mul_b'] = F(
            0x13cbd79f7463c1ca7cfd8c662c6b69d2ef3b2bf383be658c3217117a624973f3)
        test_vectors['a_div_b'] = F(
            0x3ddbe37f9605b5d0f48eabfb2bd0fbd1acd8100d6f7d1172171795192e50fc75)
        test_vectors['non_square'] = F(0x5)
        test_vectors['sqrt_b'] = F(
            0x2ca0f83ebf679a313a66801928cad698850f9b0f59a7e7f344e4929eeae32bc7)
        return F, test_vectors

    def test_random(self):
        """Check if two random elements are different.

        It happens with probability 1/`p`, small for large `p`."""
        F, test_vectors = self.set_up_field()
        a = F.random()
        b = F.random()
        self.assertFalse(a == b)

    def test_mul(self):
        """Test vector multiplication test."""
        F, test_vectors = self.set_up_field()
        a = test_vectors['a']
        b = test_vectors['b']
        a_mul_b = a*b
        self.assertEqual(a_mul_b, test_vectors['a_mul_b'])

    def test_div(self):
        """Test vector division test."""
        F, test_vectors = self.set_up_field()
        a = test_vectors['a']
        b = test_vectors['b']
        a_div_b = a/b
        self.assertEqual(a_div_b, test_vectors['a_div_b'])

    def test_sqrt(self):
        """Test vector square root test."""
        F, test_vectors = self.set_up_field()
        sq = test_vectors['b']
        root = sq.sqrt()
        self.assertEqual(root*root, sq)
        self.assertTrue(
            root == test_vectors['sqrt_b'] or root == -test_vectors['sqrt_b'])

    def test_is_square(self):
        """Legendre symbol test on small squares and a non-square."""
        F, test_vectors = self.set_up_field()
        self.assertFalse(test_vectors['non_square'].is_square())
        for i in range(F.non_square.value):
            self.assertTrue(F(i).is_square())


if __name__ == '__main__':
    unittest.main()
