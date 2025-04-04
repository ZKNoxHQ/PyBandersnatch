# -*- coding: utf-8 -*-
import unittest
from src.curve.montgomery import Montgomery
from src.field import Field
from src.primitives.xecdh import xECDH


class TestKeyExchange(unittest.TestCase):

    def set_up_key_exchange(self, secret=None):
        F = Field(
            0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001)
        a = F(0x4247698f4e32ad45a293959b4ca17afa4a2d2317e4c6ce5023e1fd63d1b5de98)
        b = F(5)
        r = 0x1cfb69d4ca675f520cce760202687600ff8f87007419047174fd06b52876e7e1
        h = 4
        E = Montgomery(a, b, r, h)
        key_exchange = xECDH(E, secret)
        return key_exchange

    def test_key_exchange(self):
        alice = self.set_up_key_exchange()
        bob = self.set_up_key_exchange()
        key1 = alice.compute_shared_secret(bob.public_key)
        key2 = bob.compute_shared_secret(alice.public_key)
        self.assertEqual(key1, key2)
