# -*- coding: utf-8 -*-
import unittest
from src.curve.montgomery import Montgomery
from src.field import Field
from src.primitives.xecdh import xECDH


class TestBandersnatchKeyExchange(unittest.TestCase):

    def set_up_key_exchange(self, secret=None, base=None):
        F = Field(
            0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001)
        a = F(0x4247698f4e32ad45a293959b4ca17afa4a2d2317e4c6ce5023e1fd63d1b5de98)
        b = F(5)
        r = 0x1cfb69d4ca675f520cce760202687600ff8f87007419047174fd06b52876e7e1
        h = 4
        E = Montgomery(a, b, r, h, g=base, glv=True)
        key_exchange = xECDH(E, secret)
        return key_exchange

    def test_key_exchange(self):
        alice = self.set_up_key_exchange()
        bob = self.set_up_key_exchange()
        key1 = alice.compute_shared_secret(bob.public_key)
        key2 = bob.compute_shared_secret(alice.public_key)
        self.assertEqual(key1, key2)


class TestEd25519KeyExchange(unittest.TestCase):

    def set_up_key_exchange(self, secret=None, base=None):
        F = Field(
            2**255-19)
        a = F(486662)
        b = F(1)
        r = 7237005577332262213973186563042994240857116359379907606001950938285454250989
        h = 8
        E = Montgomery(a, b, r, h, g=base, glv=False)
        key_exchange = xECDH(E, secret)
        return key_exchange

    def test_rfc7748__1(self):
        # Test vector 1 from RFC 7748
        u_in = "e6db6867583030db3594c1a424b15f7c726624ec26b3353b10a903a6d0ab1c4c"
        scalar = "a546e36bf0527c9d3b16154b82465edd62144c0ac1fc5a18506a2244ba449ac4"
        u_out = "c3da55379de9c6908e94ea4df28d084f32eccf03491c71f754b4075577a28552"

        g_x = int.from_bytes(bytes.fromhex(
            u_in), 'little')

        alice = self.set_up_key_exchange(
            secret=scalar,
            base=[g_x, 1]
        )
        assert alice.public_key.normalize().x == int.from_bytes(
            bytes.fromhex(u_out), 'little')

    def test_rfc7748_2(self):
        # Test vector from RFC 7748
        u_in = "e5210f12786811d3f4b7959d0538ae2c31dbe7106fc03c3efc4cd549c715a493"
        scalar = "4b66e9d4d1b4673c5ad22691957d6af5c11b6421e0ea01d42ca4169e7918ba0d"
        u_out = "95cbde9476e8907d7aade45cb4b873f88b595a68799fa152e6f8f7647aac7957"

        g_x_bytes = bytearray.fromhex(
            u_in)
        g_x_bytes[31] &= 0x7f
        g_x = int.from_bytes(g_x_bytes, 'little')
        alice = self.set_up_key_exchange(
            secret=scalar,
            base=[g_x, 1]
        )
        assert alice.public_key.normalize().x == int.from_bytes(
            bytes.fromhex(u_out), 'little')
