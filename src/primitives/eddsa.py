# -*- coding: utf-8 -*-
import hashlib
import os
from src.field import Field
from src.curve.edwards import Edwards


class EdDSA:
    def __init__(self, curve, private_key=None):
        self.curve = curve
        self.private_key = private_key or self.generate_private_key()
        self.public_key = self.generate_public_key()

    def __repr__(self):
        return "---\nEdDSA with:\n\tPrivate key: {}\n\tPublic key:{}\n---".format(self.private_key, self.public_key)

    def generate_private_key(self):
        """Generates a private key as  32 bytes (corresponds to an integer mod self.curve.r).

        WARNING: urandom is maybe not secure, but the signature is passed into sha512 later.

        """
        return os.urandom(256//8)

    def generate_public_key(self):
        """Generates the public key using the generator and the hash function sha512.

        Reference: https://datatracker.ietf.org/doc/html/rfc8032.

        """
        if len(self.private_key) != 32:
            raise Exception("Bad size of private key")
        h = hashlib.sha512(self.private_key).digest()
        a = int.from_bytes(h[:32], "little")
        a &= (1 << 254) - 8
        a |= (1 << 254)
        return (a * self.curve.generator).encode_base(256)

    def sign(self, msg):
        """Signature of a message.

        Reference: https://datatracker.ietf.org/doc/html/rfc8032.

        """
        h = hashlib.sha512(self.private_key).digest()
        a = int.from_bytes(h[:32], "little")
        a &= (1 << 254) - 8  # ensure no multiple of 4
        a |= (1 << 254)
        prefix = h[32:]
        A = (a*self.curve.generator).encode_base(256)
        r = int.from_bytes(hashlib.sha512(
            prefix + str.encode(msg)).digest(), "little") % self.curve.r
        R = r * self.curve.generator
        Rs = R.encode_base(256)
        h = int.from_bytes(hashlib.sha512(
            Rs+A+str.encode(msg)).digest(), "little") % self.curve.r
        s = (r + h * a) % self.curve.r
        return Rs + int.to_bytes(s, 32, "little")

    def verify(self, msg, signature):
        """Verification of a signature of a message.

        Reference: https://datatracker.ietf.org/doc/html/rfc8032.

        """
        if len(self.public_key) != 32:
            raise Exception("Bad public key length")
        if len(signature) != 64:
            Exception("Bad signature length")
        A = self.curve.decode_base(self.public_key, 256)
        if not A:
            return False
        Rs = signature[:32]
        R = self.curve.decode_base(Rs, 256)
        if not R:
            return False
        s = int.from_bytes(signature[32:], "little")
        if s >= self.curve.r:
            return False
        h = int.from_bytes(hashlib.sha512(
            Rs + self.public_key + str.encode(msg)).digest(), "little") % self.curve.r
        s_b_minus_h_A = self.curve.generator.multi_scalar_mul(
            s, A, self.curve.r-h)
        return s_b_minus_h_A == R  # sB-hA == R ?
