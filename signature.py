# -*- coding: utf-8 -*-
import unittest
from field import Field
from curve import Curve
import hashlib
import random
from gmpy2 import mod


class Signature:
    def __init__(self, curve, private_key=None):
        self.curve = curve
        self.private_key = private_key or self.generate_private_key()
        self.public_key = self.generate_public_key()
        self.public_key_minus_generator = self.curve.generator  # TODO

    def generate_private_key(self):
        """Generates a private key as a random integer modulo r."""
        return random.randint(1, self.curve.r - 1)

    def generate_public_key(self):
        """Generates the public key using the generator."""
        return self.private_key * self.curve.generator

    def sign(self, message):
        """Signature of a message using the private key."""
        e = int(hashlib.sha256(message.encode()).hexdigest(), 16)  # Hash message
        while True:
            k = random.randint(1, self.curve.r - 1)
            p = k * self.curve.generator
            r = p.x/p.z
            if r == 0:
                continue

            s = (e + r * self.private_key)/k
            if s != 0:
                break
        # TODO is the modulo here or somewhere else?
        return (mod(r.value, self.curve.r), mod(s.value, self.curve.r))

    def verify(self, message, signature):
        """Verification of a signature for a given message."""
        r, s = signature

        if not (1 < r and r < self.curve.r and 1 < s and s < self.curve.r):
            return False

        e = int(hashlib.sha256(message.encode()).hexdigest(), 16)

        Fr = Field(
            0x1cfb69d4ca675f520cce760202687600ff8f87007419047174fd06b52876e7e1)

        w = Fr(1)/Fr(s)
        u1 = Fr(e)*w
        u2 = Fr(r)*w

        p = self.curve.generator.multi_scalar_mul(
            u1.value, self.public_key, u2.value, self.public_key_minus_generator)
        return (p.x/p.z).value == r


class TestSignature(unittest.TestCase):

    def set_up_signature(self, secret=None):
        F = Field(
            0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001)
        a = F(0x4247698f4e32ad45a293959b4ca17afa4a2d2317e4c6ce5023e1fd63d1b5de98)
        b = F(5)
        r = 0x1cfb69d4ca675f520cce760202687600ff8f87007419047174fd06b52876e7e1
        h = 4
        E = Curve(a, b, r, h)
        signature = Signature(E, secret)
        return signature

    def test_to_start(self):
        signature = self.set_up_signature()
        signature = self.set_up_signature(123)
        self.assertEqual(1, 1)

    def test_signature_generation(self):
        signature = self.set_up_signature()
        print("pk = {}".format(signature.public_key))
        print("sk = {}".format(signature.private_key))
        (r, s) = signature.sign("Salut")
        print("signature = ({}, {})".format(r, s))
        print(signature.verify("Salut", (r, s)))


if __name__ == '__main__':
    unittest.main()
