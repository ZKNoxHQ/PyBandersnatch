# -*- coding: utf-8 -*-
import unittest
from curve import Curve
from field import Field
from key_exchange import KeyExchange


class TestKeyExchange(unittest.TestCase):

    def set_up_key_exchange(self, secret=None):
        F = Field(
            0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001)
        a = F(0x4247698f4e32ad45a293959b4ca17afa4a2d2317e4c6ce5023e1fd63d1b5de98)
        b = F(5)
        r = 0x1cfb69d4ca675f520cce760202687600ff8f87007419047174fd06b52876e7e1
        h = 4
        E = Curve(a, b, r, h)
        key_exchange = KeyExchange(E, secret)
        return key_exchange

    def test_key_exchange(self):
        alice = self.set_up_key_exchange()
        bob = self.set_up_key_exchange()
        key1 = alice.compute_shared_secret(bob.public_key)
        key2 = bob.compute_shared_secret(alice.public_key)
        self.assertEqual(key1, key2)

    def run_all_test(self):
        print("Tests for KeyExchange")
        print("------")
        for method_name in dir(self):
            if method_name.startswith("test_"):
                method = getattr(self, method_name)
                if callable(method):
                    print(f"Running: {method_name}")
                    method()
        print("------")
        print()
