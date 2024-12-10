# -*- coding: utf-8 -*-
import random
from src.field import Field
from src.curve.montgomery import Montgomery


class xECDH:
    def __init__(self, curve, private_key=None):
        self.curve = curve
        self.private_key = private_key or self.generate_private_key()
        self.public_key = self.generate_public_key()

    def generate_private_key(self):
        """Generates a private key as a random integer modulo r.

        WARNING: not secure.

        """
        return random.randint(1, self.curve.r - 1)

    def generate_public_key(self):
        """Generates the public key using the generator."""
        return self.private_key * self.curve.g

    def compute_shared_secret(self, other_public_key):
        """Compute the shared secret `secret_key` * `other_public_key`."""
        shared_secret = self.private_key * other_public_key
        return shared_secret
