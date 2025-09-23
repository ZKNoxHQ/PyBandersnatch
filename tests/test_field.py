# tests/test_field.py
import unittest
import os


class FieldBase:
    VECTOR_FILE = None

    def set_up_field(self):
        path = os.path.join("tests/vectors", self.VECTOR_FILE)
        try:
            with open(path, "r") as file:
                exec(file.read(), globals())
        except FileNotFoundError:
            raise unittest.SkipTest(
                f"File '{path}' not found. Please generate it with Sage."
            )
        return F, test_vectors  # type: ignore

    def test_random(self):
        F, _ = self.set_up_field()
        a, b = F.random(), F.random()
        self.assertNotEqual(a, b)

    def test_mul(self):
        F, tv = self.set_up_field()
        self.assertEqual(tv['a'] * tv['b'], tv['a_mul_b'])

    def test_div(self):
        F, tv = self.set_up_field()
        self.assertEqual(tv['a'] / tv['b'], tv['a_div_b'])

    def test_sqrt(self):
        F, tv = self.set_up_field()
        sq = tv['b']
        root = sq.sqrt()
        self.assertEqual(root * root, sq)
        self.assertTrue(root == tv['sqrt_b'] or root == -tv['sqrt_b'])

    def test_is_square(self):
        F, tv = self.set_up_field()
        self.assertFalse(tv['non_square'].is_square())
        for i in range(F.non_square.value):
            self.assertTrue(F(i).is_square())


class TestBandersnatchField(FieldBase, unittest.TestCase):
    VECTOR_FILE = "bandersnatch_field.py"


class TestEd25519Field(FieldBase, unittest.TestCase):
    VECTOR_FILE = "ed25519_field.py"
