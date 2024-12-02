# -*- coding: utf-8 -*-
import unittest
from src.curve.edwards import Edwards
from src.field import Field
from src.primitives.eddsa import EdDSA


class TestEdDSA(unittest.TestCase):

    def set_up_eddsa(self, secret=None):
        F = Field(
            0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001)
        a = F(-5)
        d = F(45022363124591815672509500913686876175488063829319466900776701791074614335719)
        r = 0x1cfb69d4ca675f520cce760202687600ff8f87007419047174fd06b52876e7e1
        h = 4
        E = Edwards(a, d, r, h)
        eddsa = EdDSA(E, private_key=secret)
        return eddsa

    def test_sign_verify(self):
        """Signature verification works"""
        alice = self.set_up_eddsa(secret=b"Que boludo... my llave es fija!!")
        sig_1 = alice.sign("Buenas che")
        assert alice.verify("Buenas che", sig_1)

        bob = self.set_up_eddsa()
        sig_2 = bob.sign("Boa noite")
        assert bob.verify("Boa noite", sig_2)
