# -*- coding: utf-8 -*-
import unittest
from src.curve.edwards import Edwards
from src.field import Field
from src.primitives.eddsa import EdDSA
from tests.test_edwards_25519 import TestEdwards25519


class TestEdDSA25519(unittest.TestCase):

    def set_up_eddsa(self, secret=None):
        E, _ = TestEdwards25519.set_up_curve(TestEdwards25519)
        return EdDSA(E, private_key=secret)

    def test_sign_verify(self):
        """Signature verification works"""
        alice = self.set_up_eddsa(secret=b"Que boludo... my llave es fija!!")
        sig_1 = alice.sign("Buenas che")
        assert alice.verify("Buenas che", sig_1)

        bob = self.set_up_eddsa()
        sig_2 = bob.sign("Boa noite")
        assert bob.verify("Boa noite", sig_2)
