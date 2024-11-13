import random
import unittest
from field import Field
from curve import Curve


class KeyExchange:
    def __init__(self, curve, private_key=None):
        self.curve = curve
        self.private_key = private_key or self.generate_private_key()
        self.public_key = self.generate_public_key()

    def generate_private_key(self):
        """Generates a private key as a random integer modulo r."""
        return random.randint(1, self.curve.r - 1)

    def generate_public_key(self):
        """Generates the public key using the generator."""
        return self.private_key * self.curve.generator

    def compute_shared_secret(self, other_public_key):
        """Compute the shared secret `secret_key` * `other_public_key`."""
        shared_secret = self.private_key * other_public_key
        return shared_secret


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


if __name__ == '__main__':
    unittest.main()
