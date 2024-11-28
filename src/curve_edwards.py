# -*- coding: utf-8 -*-
import argparse
from field import Field
from gmpy2 import mpq, mpz, mod


class CurveEdwards:
    def __init__(self, a, d, r, h):
        self.field = a.field
        self.a = a
        self.d = d
        self.r = r
        self.h = h
        self.a24 = (self.a+2)/4
        self.generator = self.Point(self.field(1), self.field(
            0x39c8932646733e6987574d2c518adfb172294f2c8c95981c0c95f8188762dd7a), self.field(1), self)

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
            return self.z == 0 or other.z == 0 or (self.x * other.z == other.x * self.z and self.y * other.z == other.y * self.z)

        def normalize(self):
            """Affine representation of the projective point."""
            if self.z == 0:
                return self.curve(0, 1, 0)
            return self.curve(self.x/self.z, self.y/self.z, 1)

        def neg(self):
            return self.curve(-self.x, self.y, self.z)

        def in_curve(self):
            """Returns the curve membership boolean."""
            if self.z == 0:
                return True
            x = self.x/self.z
            a = self.curve.a
            d = self.curve.d
            return ((1-a*x**2)/(1-d*x**2)).is_square()

        def dbl(self):
            """Doubling algorithm.

            Reference:
            https://eprint.iacr.org/2008/013.pdf page 12.

            """
            if self.z == 0:
                return self.curve(0, 1, 0)
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
            if self.z == 0:
                return q
            if q.z == 0:
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
                return self.curve(0, 1, 0)
            if k < 0:
                k = -k
                self = self.neg()
            res = self.curve(0, 1, 0)
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
            if res.z == 0:
                return self.curve(1, 0, 0)
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
                return self.curve(0, 1, 0)

            if s1 > s0:
                s0, p0, s1, p1 = s1, p1, s0, p0
            # s1 ≤ s0

            res = self.curve(0, 1, 0)
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

        def is_prime_order(self, N):
            """Returns the boolean corresponding to `self.order() == N`."""
            return self.naive_mul(N).x == 0

        def φ(self):
            """Endomorphism sqrt(-2).

            TODO can be optimized with factorization of multi-variable polynomials.
            Obtained from `sage tests/φ_edwards.sage`.

            """
            ay4 = 0x5d149954d89e9ae0796a0a7bcd57a1c5a8cb126b81e2333f8b1c3362af1673f3
            ay2z2 = 0x471115903c695fcf3153011b2354fe7fc1f7385cb321eaf52e422633906c1199
            az4 = 0x36eef32fe7d73d4462bbecd2e886a5b521811d030558648d889f592edf5fb3d1
            by2 = 0x4b37256660a3d1343b42833f3d7e7cc77921a3766303c1b1d00a4126765bff2e
            bz2 = 0x721f22321a0a7afe7c81a39b08e55cd926a687b17dde811c7ddbb0f949dca6c0
            cy2 = 0x20b21e58881722d68c92fa09709ea65d716e869843e94c821df033483694a51a
            cz2 = 0x28b681ecc8f9ac13f7f754c8cc235b3dda9c008c9cfa9a4d2ff5bed889a400d3

            x = self.x
            y = self.y
            z = self.z
            x_r = x * (ay4 * y**4 + ay2z2 * y**2 * z**2 + az4 * z**4)
            y_r = y * z**2 * (by2 * y**2 + bz2 * z**2)
            z_r = y * z**2 * (cy2 * y**2 + cz2 * z**2)
            return self.curve(x_r, y_r, z_r)

        def glv(self, k, constant_time=False):
            """GLV scalar multiplication `k`*`self`.

            A constant time option is available.
            Reference:
            https://www.iacr.org/archive/crypto2001/21390189.pdf
            Obtained from `sage tests/φ_edwards.sage`.

            """
            if k == 0:
                return self.curve(0, 1, 0)

            M1 = [-113482231691339203864511368254957623327,
                  10741319382058138887739339959866629956]
            M2 = [21482638764116277775478679919733259912,
                  113482231691339203864511368254957623327]
            N1 = [-113482231691339203864511368254957623327,
                  10741319382058138887739339959866629956]

            b = [round(mpq(k*N1[0], self.curve.r)),
                 round(mpq(k*N1[1], self.curve.r))]
            k1 = k-b[0] * M1[0] - b[1] * M2[0]
            k2 = -b[0] * M1[1] - b[1] * M2[1]
            return self.multi_scalar_mul(k1, self.φ(), k2)

        def __rmul__(self, k):
            """Scalar multiplication with the scalar give first.

            Computed using GLV.

            """
            return self.glv(k)
