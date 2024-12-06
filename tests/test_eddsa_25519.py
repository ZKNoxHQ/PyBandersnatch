# -*- coding: utf-8 -*-
import unittest
from src.curve.edwards25519 import Edwards25519
from src.field import Field
from src.primitives.eddsa import EdDSA


class TestEdDSA(unittest.TestCase):

    def set_up_eddsa(self, secret=None):
        F = Field(
            0x7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffed)
        a = F(57896044618658097711785492504343953926634992332820282019728792003956564819944)
        d = F(11790395817372903580333940030740964167805592400755249022757551653560006957928)
        r = 7237005577332262213973186563042994240857116359379907606001950938285454250989
        h = 8
        E = Edwards25519(a, d, r, h)
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
