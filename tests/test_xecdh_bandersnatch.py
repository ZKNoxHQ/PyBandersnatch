# -*- coding: utf-8 -*-
import unittest
from src.curve.montgomery import Montgomery
from src.field import Field
from src.primitives.xecdh import xECDH
from tests.test_montgomery_bandersnatch import TestMontgomeryBandersnatch


class TestxECDHBandersnatch(unittest.TestCase):

    def set_up_key_exchange(self, secret=None):
        E, _ = TestMontgomeryBandersnatch.set_up_curve(
            TestMontgomeryBandersnatch)
        return xECDH(E, secret)

    def test_key_exchange(self):
        alice = self.set_up_key_exchange()
        bob = self.set_up_key_exchange()
        key1 = alice.compute_shared_secret(bob.public_key)
        key2 = bob.compute_shared_secret(alice.public_key)
        self.assertEqual(key1, key2)
