# -*- coding: utf-8 -*-
from gmpy2 import random_state, mpz_random, invert, mpz, powmod, sign, mod


class Field():

    def __init__(self, p):
        self.p = mpz(p)
        self._rand_state = random_state()

        # 2-adicity
        self.two_adicity = 0
        while (p-1) % 2**self.two_adicity == 0:
            self.two_adicity += 1
        self.two_adicity -= 1

        # Quadratic non-residue
        self.non_square = self.Element(mpz(1), self)
        while self.non_square.is_square():
            self.non_square += self.Element(1, self)

    def __str__(self):
        return "Finite field of characteristic {}".format(self.p)

    def __call__(self, value):
        return self.Element(value, self)

    def random(self):
        """Compute a random element of `self`."""
        # probably not secure
        return self.Element(mpz_random(self._rand_state, self.p), self)

    class Element:
        def __init__(self, value, field):
            self.value = value % field.p
            self.field = field

        def __eq__(self, other):
            """Return the equality boolean between `self` and `other`."""
            if isinstance(other, int):
                return self.value == other % self.field.p
            elif isinstance(other, self.field.Element) and self.field.p == other.field.p:
                return self.value == other.value
            raise TypeError("Cannot add elements from different fields.")

        def __add__(self, other):
            """Addition of `self` and `other`."""
            other = self.field(other) if isinstance(
                other, int) or isinstance(other, mpz) else other
            if self.field.p == other.field.p:
                return self.field(self.value + other.value)
            raise TypeError("Cannot add elements from different fields.")

        def __radd__(self, other):
            """Addition when `other` is given first (mostly for `int` type)."""
            return self+other

        def __neg__(self):
            """Negation of `self`."""
            return self.field(-self.value % self.field.p)

        def __sub__(self, other):
            """Difference of `self` and `other`."""
            other = self.field(other) if isinstance(
                other, int) or isinstance(other, mpz) else other
            if self.field.p == other.field.p:
                return self.field(self.value - other.value)
            raise TypeError("Cannot subtract elements from different fields.")

        def __rsub__(self, other):
            """Substraction when `other` is given first (mostly for `int` type)."""
            return -self+other

        def __mul__(self, other):
            """Multiplication of `self` and `other`."""
            other = self.field(other) if isinstance(
                other, int) or isinstance(other, mpz) else other
            if self.field.p == other.field.p:
                return self.field((self.value * other.value) % self.field.p)
            raise TypeError("Cannot multiply elements from different fields.")

        def __rmul__(self, other):
            """Multiplication when `other` is given first (mostly for `int` type)."""
            return self * other

        def __pow__(self, exponent):
            """Modular exponentiation `self` to the power `exponent`."""
            result = powmod(self.value, exponent, self.field.p)
            return self.field(result)

        def is_square(self):
            """Legendre symbol of `self`."""
            # TODO optimize the case (p,q) = (q,p) * something
            if self.value == 0:
                return True
            else:
                return self ** ((self.field.p-1) >> 1) == 1

        def __truediv__(self, other):
            """Division of `self` by `other` modulo `self.field.p`."""
            if isinstance(other, int):
                return self/self.field(other)
            if isinstance(other, self.field.Element) and self.field == other.field:
                inverse_other = invert(other.value, self.field.p)
                return self.field((self.value * inverse_other) % self.field.p)
            # TODO division by zero error?
            raise TypeError("Cannot divide elements from different fields.")

        def __repr__(self):
            return f"{self.value}"

        def sqrt(self):
            """Square root of `self`.

            Computed using Tonelli-Shanks algorithm.
            Assumes that `self` is a square.
            Details here:
            https://www.lvzl.fr/teaching/2020-21/docs/AA-devoir-1-sujet.pdf
            https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm

            """

            if self.value == 0:
                return self

            if self.field.p % 4 == 3:
                # TODO ASSERT IS SQUARE BEFORE RETURN
                return powmod(self.value, (self.field.p+1)//4, self.field.p)

            # Tonelli-Shanks
            z = self.field.non_square
            M = self.field.two_adicity
            Q = (self.field.p-1)//2**M

            c = z**Q
            t = self**Q
            r = self**((Q+1) >> 1)
            # TODO REMOVE THIS CHECK
            assert (t**(1 << self.field.two_adicity) == 1)

            while t != 1:
                # Find the least i (0 < i < m) such that t^(2^i) ≡ 1 (mod p)
                t2i = t
                i = 0
                while t2i != 1:
                    t2i = t2i**2
                    i += 1

                # Update variables for the next iteration
                b = c**(1 << (M-i-1))
                M = i
                c = b*b
                t = t*c
                r = r*b

            return r if t == 1 else None
