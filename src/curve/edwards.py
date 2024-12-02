# -*- coding: utf-8 -*-
from math import floor
from src.field import Field
from gmpy2 import mpq, mpz, mod


class Edwards:
    def __init__(self, a, d, r, h):
        self.field = a.field
        self.a = a
        self.d = d
        self.r = r
        self.h = h
        self.generator = self.Point(self.field(3), self.field(
            0x2d418cc584d9c9df8750a436fac98068949d14c7bdce4034fe792e4c14e30a3f), self.field(1), self)

    def __repr__(self):
        return "Edwards curve defined by {}*x^2 + y^2 = 1 + {} * x^2*y^2".format(self.a, self.d)

    __str__ = __repr__

    def __call__(self, x, y, z):
        if isinstance(x, int) or isinstance(x, mpz):
            x = self.field(x)
        if isinstance(y, int) or isinstance(y, mpz):
            y = self.field(y)
        if isinstance(z, int) or isinstance(z, mpz):
            z = self.field(z)
        return self.Point(x, y, z, self)

    def random(self):
        """Returns a random point of `self`."""
        x = self.field.random()
        while not ((1-self.a*x**2)/(1-self.d*x**2)).is_square():
            x = self.field.random()
        y = ((1-self.a*x**2)/(1-self.d*x**2)).sqrt()
        return self.Point(x, y, self.field(1), self)

    def j_inv(self):
        """Returns the j-invariant of `self`.

        Reference:
        https://eprint.iacr.org/2008/013.pdf page 3.

        """
        a = self.a
        d = self.d
        return 16*(a**2 + 14*a*d + d**2)**3 / (a*d*(a-d)**4)

    def decode_base(self, s, b):
        """Decoding following the formate of RFC 8032.

        Reference: https://datatracker.ietf.org/doc/html/rfc8032

        """
        # Check that point encoding is the correct length.
        if len(s) != b//8:
            return (None, None)
        # Extract signbit.
        xs = s[(b-1)//8] >> ((b-1) & 7)
        # Decode y.  If this fails, fail.
        p = self.field.p
        rv = int.from_bytes(s, byteorder="little") % (2**(b-1))
        y = self.field(rv) if rv < p else None
        if y is None:
            return (None, None)
        # ax² + y² = 1 + dx²y² => x² = (1-y²)/(a -dy²)
        # Try to recover x.
        # If it does not exist, or if zero and xs are wrong, fail.
        x = ((1-y**2)/(self.a-self.d*y**2)).sqrt()
        if x is None or (x == 0 and xs != (x.value % p) % 2):
            return (None, None)
        # If sign of x isn't correct, flip it.
        if (x.value % p) % 2 != xs:
            x = -x
        # Return the constructed point.
        return self(x, y, 1)

    class Point:
        def __init__(self, x, y, z, curve):
            self.x = x
            self.y = y
            self.z = z
            self.curve = curve

        def __repr__(self):
            return "Point ({}, {}, {})".format(self.x, self.y, self.z)

        def __eq__(self, other):
            """Return the equality boolean between `self` and `other`.

            Points are represented in projective coordinates.
            """
            if self.z != 0 and other.z != 0:
                return (self.x * other.z == other.x * self.z and self.y * other.z == other.y * self.z)
            else:
                if self.z == 0 and other.z == 0:
                    # self and other are both at infinity
                    if self.x == 0 and other.x == 0:
                        # (0,a,0) == (0,b,0)
                        return True
                    elif self.y == 0 and other.y == 0:
                        # (a,0,0) == (b,0,0)
                        return True
                    else:
                        # (a,0,0) ≠ (0,b,0)
                        return False
                    return self.x
                else:
                    # only one of the two points is at infinity
                    return False

        def normalize(self):
            """Affine representation of the projective point."""
            if self.z == 0:
                if self.x == 0:
                    return self.curve(0, 1, 0)
                elif self.y == 0:
                    return self.curve(1, 0, 0)
                raise ("This should not happen")
            return self.curve(self.x/self.z, self.y/self.z, 1)

        def neg(self):
            return self.curve(-self.x, self.y, self.z)

        def in_curve(self):
            """Returns the curve membership boolean."""
            if self.z == 0:
                return self.x*self.y == 0
            x = self.x/self.z
            a = self.curve.a
            d = self.curve.d
            return ((1-a*x**2)/(1-d*x**2)).is_square()

        def dbl(self):
            """Doubling algorithm.

            Reference:
            https://eprint.iacr.org/2008/013.pdf page 12.

            """
            if self.z == 0 and self.x * self.y == 0:  # (1,0,0) and (0,1,0) are of order 2
                return self.curve(0, 1, 1)
            if self == self.curve(0, 1, 1):
                return self
            x = self.x
            y = self.y
            z = self.z
            b = (x+y)**2
            c = x**2
            d = y**2
            e = self.curve.a * c
            f = e + d
            h = z**2
            j = f-2*h
            x_r = (b-c-d)*j
            y_r = f * (e-d)
            z_r = f*j
            return self.curve(x_r, y_r, z_r)

        def add(self, q):
            """Addition algorithm.

            Reference:
            https://eprint.iacr.org/2008/013.pdf page 12.

            """
            if self == self.curve(0, 1, 1):
                return q
            if q == self.curve(0, 1, 1):
                return self
            x_p, y_p, z_p = self.x, self.y, self.z
            x_q, y_q, z_q = q.x, q.y, q.z
            a = z_p * z_q
            b = a**2
            c = x_p * x_q
            d = y_p * y_q
            e = self.curve.d * c * d
            f = b-e
            g = b+e
            x_r = a*f*((x_p+y_p) * (x_q+y_q) - c - d)
            y_r = a*g*(d-self.curve.a*c)
            z_r = f*g
            return self.curve(x_r, y_r, z_r)

        def naive_mul(self, k):
            """Scalar multiplication `k` * `self`.

            Double-and-add algorithm
            From most significant bit (MSB) to least significant bit (LSB)
            TODO not constant time

            """
            if k == 0:
                return self.curve(0, 1, 1)
            if k < 0:
                k = -k
                self = self.neg()
            res = self.curve(0, 1, 1)
            temp = self
            k = int(k)
            n = k.bit_length()
            # msb -> lsb
            while n > 0:
                res = res.dbl()
                n -= 1
                if (k >> n) & 1:
                    res = res.add(self)
                k &= ~(1 << n)
            return res

        def multi_scalar_mul(self, k1, other, k2):
            """Multi scalar multiplication `k1` * `self` + `k2` * `other`.

            From most significant bit (MSB) to least significant bit (LSB)
            TODO not constant time.

            """
            s0, s1, p0, p1 = k1, k2, self, other

            p0 = p0.neg() if int(s0) < 0 else p0
            p1 = p1.neg() if int(s1) < 0 else p1
            s0 = abs(int(s0))
            s1 = abs(int(s1))

            if s0 == 0 and s1 == 0:
                return self.curve(0, 1, 1)

            if s1 > s0:
                s0, p0, s1, p1 = s1, p1, s0, p0
            # s1 ≤ s0

            res = self.curve(0, 1, 1)
            prec = [[res, p1], [p0, p0.add(p1)]]
            s0 = int(s0)
            s1 = int(s1)
            n = s0.bit_length()
            while n > 0:
                res = res.dbl()
                n -= 1
                res = res.add(prec[(s0 >> n) & 1][(s1 >> n) & 1])
                s0 &= ~(1 << n)
                s1 &= ~(1 << n)
            return res

        def is_prime_order(self, n):
            """Returns the boolean corresponding to `self.order() == n`.

            Bandersnatch has only points of order 2 or r (or infinity).
            For other curves, other cases could happen.

            """
            if self == self.curve(0, 1, 1):
                return n == 0
            n_times_p = self.naive_mul(n)
            if n == 2:  # TODO: does it work in the general case?
                return n_times_p.z == 0 or self == self.curve(0, -1, 1)
            else:
                return n_times_p == self.curve(0, 1, 1)

        def φ(self):
            """Endomorphism sqrt(-2).

            TODO can be optimized with factorization of multi-variable polynomials.
            Obtained using `sage sage/φ.sage`.

            """
            ay4 = 0x1d46e71b2d28e06c42bc1f5a41f4a0156d070863689e8862eb12927f72f308c3
            ay2z2 = 0x20b21e58881722d68c92fa09709ea65d716e869843e94c821df033483694a51a
            az4 = 0x1373fe65dcb354e5209f902de5b37008d6c2721d8d6d5fb556e8b7e969c053c9
            by2 = 0x33937d60e9a0dd55ed1f9030e7c8b6fa9c42e1e41d2f1361a0fed9630f711cae
            bz2 = 0x39a33e54438fe0155ae18e93205d4395acfe4be1127ca6458fc0270450b1b50d
            cy2 = 0x2cdc91c2ed341d7901e6d6ece64cd98591c66ba64cdc7109d1bdd9cb6f93ee68
            cz2 = 0x405a29f23ffc9ff2461a47d721d9210ab77ac21ee2cf489d5f01269bf08ee353

            x = self.x
            y = self.y
            z = self.z
            x_r = x * (ay4 * y**4 + ay2z2 * y**2 * z**2 + az4 * z**4)
            y_r = y * z**2 * (by2 * y**2 + bz2 * z**2)
            z_r = y * z**2 * (cy2 * y**2 + cz2 * z**2)
            return self.curve(x_r, y_r, z_r)

        def glv(self, k):
            """GLV scalar multiplication `k`*`self`.

            WARNING: this does not work when k = r for example!!!!!!!
            Reference:
            https://www.iacr.org/archive/crypto2001/21390189.pdf

            """
            if k == 0:
                return self.curve(0, 1, 1)

            m1 = int(-113482231691339203864511368254957623327)
            m2 = int(10741319382058138887739339959866629956)
            m3 = int(21482638764116277775478679919733259912)
            b = [floor(mpq(k*m1, self.curve.r)),
                 floor(mpq(k*m2, self.curve.r))]
            k1 = k-b[0] * m1 - b[1] * m3
            k2 = -b[0] * m2 - b[1] * -m1
            return self.multi_scalar_mul(k1, self.φ(), k2)

        def __rmul__(self, k):
            """Scalar multiplication with the scalar give first.

            Computed using GLV.

            """
            return self.glv(k)

        def encode_base(self, b):
            """Decoding following the formate of RFC 8032.

            Reference: https://datatracker.ietf.org/doc/html/rfc8032

            """
            xp, yp = self.x/self.z, self.y/self.z
            p = self.curve.field.p
            s = bytearray(int(yp.value % p).to_bytes(b//8, byteorder='little'))
            if (xp.value % p) % 2 != 0:
                s[(b-1)//8] |= 1 << (b-1) % 8
            return s
