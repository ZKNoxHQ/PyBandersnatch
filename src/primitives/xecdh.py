# -*- coding: utf-8 -*-
import random
from src.field import Field
from src.curve.montgomery import Montgomery


def decode_scalar_x25519(hex_string: str) -> int:
    """Decode X25519 scalar according to RFC 7748.
    """
    k_bytes = bytearray.fromhex(hex_string)
    k_bytes[0] &= 248   # Clear bits 0,1,2
    k_bytes[31] &= 127  # Clear bit 255
    k_bytes[31] |= 64   # Set bit 254
    # Convert bytes to integer (little-endian)
    return int.from_bytes(k_bytes, 'little')


class xECDH:
    def __init__(self, curve, private_key=None):
        self.curve = curve
        if private_key == None:
            self.private_key = self.generate_private_key()
        else:
            self.private_key = decode_scalar_x25519(private_key)
        self.public_key = self.generate_public_key()

    def generate_private_key(self):
        """Generates a private key as a random integer modulo r.

        WARNING: not secure.

        """
        return random.randint(1, self.curve.r - 1)

    def generate_public_key(self):
        """Generates the public key using the generator."""
        return self.private_key * self.curve.generator

    def compute_shared_secret(self, other_public_key):
        """Compute the shared secret `secret_key` * `other_public_key`."""
        shared_secret = self.private_key * other_public_key
        return shared_secret
