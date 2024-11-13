# -*- coding: utf-8 -*-
import unittest
from field import Field
from curve import Curve
import hashlib
import random
from gmpy2 import mod, invert


class Signature:
    def __init__(self, curve, private_key=None):
        self.curve = curve
        self.private_key = private_key or self.generate_private_key()
        self.public_key = self.generate_public_key()
        # self.public_key_minus_generator = self.curve.generator.slow_add(
        #     self.public_key)

    def generate_private_key(self):
        """Generates a private key as a random integer modulo r."""
        return random.randint(1, self.curve.r - 1)

    def generate_public_key(self):
        """Generates the public key using the generator."""
        return self.private_key * self.curve.generator

    def sign(self, message):
        """Signature of a message using the private key."""
        e = int(hashlib.sha256(message.encode()).hexdigest(), 16)  # Hash message
        e = mod(e, self.curve.r)
        while True:
            k = 7281454368565165140637308125972207708798847059953356515906353705248896425334
            # k = random.randint(1, self.curve.r - 1)
            p = k * self.curve.generator
            print(p)
            r = (p.x/p.z).value
            if r == 0:
                continue

            s = mod((e + mod(r, self.curve.r) * self.private_key) *
                    invert(k, self.curve.r), self.curve.r)
            if s != 0:
                break

        return (r, s)  # r is not in the scalar field!

    # def verify(self, message, signature):
    #     """Verification of a signature for a given message."""

    #     assert (self.public_key.is_prime_order(self.curve.r))

    #     r, s = signature

    #     if not (1 < r and r < self.curve.r and 1 < s and s < self.curve.r):
    #         return False

    #     e = int(hashlib.sha256(message.encode()).hexdigest(), 16)
    #     e = mod(e, self.curve.r)

    #     w = invert(s, self.curve.r)  # = k/(e+r*sk)
    #     u1 = mod(e * w, self.curve.r)  # = e k / (e+r*sk)
    #     u2 = mod(r*w, self.curve.r)  # = r k / (e+r*sk)

    #     print("u1={}".format(u1))
    #     print("u2={}".format(u2))
    #     print("r={}".format(self.curve.r))

    #     p_check = (u1 * self.curve.generator).slow_add(u2*self.public_key)
    #     p = self.curve.generator.multi_scalar_mul(
    #         # (1/(e+r*sk)) (e + r sk) kG
    #         u1, self.public_key, u2, self.public_key_minus_generator)
    #     print("p1 = {}".format(p_check))
    #     print("p2 = {}:".format(p.normalize()))
    #     print("r = {}".format(r))
    #     return mod((p.x/p.z).value, self.curve.field.p) == r


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

    def test_signature_generation(self):
        signature = self.set_up_signature()
        (r, s) = signature.sign("Salut")
        f = open("example.sage", "w")
        f.write("public_key = {}\n".format(signature.public_key.normalize().x))
        f.write("signature = ({}, {})\n".format(r, s))
        f.write("message = \"{}\"".format("Salut"))
        f.close()


if __name__ == '__main__':
    unittest.main()
