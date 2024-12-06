# -*- coding: utf-8 -*-
from math import floor
from src.field import Field
from gmpy2 import mpq, mpz, mod


class Edwards:
    def __init__(self, a, d, r, h, φ_coefficients, glv_params):
        self.field = a.field
        self.a = a
        self.d = d
        self.r = r
        self.h = h
        self.g = self.generator()
        self.φ_coefficients = φ_coefficients
        self.glv_params = glv_params

    def generator(self):
        """Generator with small x-coordinate."""
        x = 0
        g = self(self.field(0), self.field(1), self.field(1))
        while not (g.is_prime_order(self.r)):
            x += 1
            while not ((1-self.a*x**2)/(1-self.d*x**2)).is_square():
                x = -x
                if x > 0:
                    x += 1
            y = ((1-self.a*x**2)/(1-self.d*x**2)).sqrt()
            # TODO: what about -y?
            g = self(self.field(x), self.field(y.value), self.field(1))
        return g

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

        def __add__(self, q):
            return self.add(q)

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
                    res += self
                k &= ~(1 << n)
            return res

        def multi_scalar_mul_2(self, k, q, l):
            """Compute k*self + l*q using GLV trick. It is a four-dimensional scalar multiplication."""
            if k == 0 and l == 0:
                return self.curve(0, 1, 1)

            m1, m2, m3 = self.curve.glv_params
            b = [floor(mpq(k*m1, self.curve.r)),
                 floor(mpq(k*m2, self.curve.r))]
            k1 = k-b[0] * m1 - b[1] * m3
            k2 = -b[0] * m2 - b[1] * -m1

            b = [floor(mpq(l*m1, self.curve.r)),
                 floor(mpq(l*m2, self.curve.r))]
            l1 = l-b[0] * m1 - b[1] * m3
            l2 = -b[0] * m2 - b[1] * -m1

            return self.multi_scalar_mul_4(k1, self.φ(), k2, q, l1, q.φ(), l2)

        def multi_scalar_mul_4(self, k1, q, k2, r, k3, s, k4):
            """Multi scalar multiplication `k1` * `self` + `k2` * `q` + `k3` * `r` + `k4` * `s`.

            From most significant bit (MSB) to least significant bit (LSB)
            TODO not constant time.

            """
            s0, s1, s2, s3, p0, p1, p2, p3 = k1, k2, k3, k4, self, q, r, s

            p0 = p0.neg() if int(s0) < 0 else p0
            p1 = p1.neg() if int(s1) < 0 else p1
            p2 = p2.neg() if int(s2) < 0 else p2
            p3 = p3.neg() if int(s3) < 0 else p3
            s0 = abs(int(s0))
            s1 = abs(int(s1))
            s2 = abs(int(s2))
            s3 = abs(int(s3))

            if s0 == 0 and s1 == 0 and s2 == 0 and s3 == 0:
                return self.curve(0, 1, 1)

            if s1 > s0:
                s0, p0, s1, p1 = s1, p1, s0, p0
            if s2 > s0:
                s0, p0, s2, p2 = s2, p2, s0, p0
            if s3 > s0:
                s0, p0, s3, p3 = s3, p3, s0, p0
            # s1,s2,s3 ≤ s0

            res = self.curve(0, 1, 1)
            prec = [[[[None]*2 for _ in range(2)]
                     for _ in range(2)] for _ in range(2)]
            prec[0][0][0][0] = res
            prec[0][0][0][1] = p3
            prec[0][0][1][0] = p2
            prec[0][0][1][1] = p2 + p3
            prec[0][1][0][0] = p1
            prec[0][1][0][1] = p1 + p3
            prec[0][1][1][0] = p1 + p2
            prec[0][1][1][1] = prec[0][0][1][1] + p1
            prec[1][0][0][0] = p0
            prec[1][0][0][1] = p0 + p3
            prec[1][0][1][0] = p0 + p2
            prec[1][0][1][1] = prec[0][0][1][1] + p0
            prec[1][1][0][0] = p0 + p1
            prec[1][1][0][1] = prec[1][0][0][1] + p1
            prec[1][1][1][0] = prec[1][1][0][0] + p2
            prec[1][1][1][1] = prec[1][1][1][0] + p3
            s0 = int(s0)
            s1 = int(s1)
            s2 = int(s2)
            s3 = int(s3)
            n = s0.bit_length()
            while n > 0:
                res = res.dbl()
                n -= 1
                res += prec[(s0 >> n) & 1][(s1 >> n) &
                                           1][(s2 >> n) & 1][(s3 >> n) & 1]
                s0 &= ~(1 << n)
                s1 &= ~(1 << n)
                s2 &= ~(1 << n)
                s3 &= ~(1 << n)
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
            """Endomorphism sqrt(-2) (or Identity if large discriminant).

            TODO can be optimized with factorization of multi-variable polynomials.
            Obtained using `sage sage/φ.sage`.

            """
            ay4, ay2z2, az4, by2, bz2, cy2, cz2, cyz = self.curve.φ_coefficients
            x = self.x
            y = self.y
            z = self.z
            x_r = x * (ay4 * y**4 + ay2z2 * y**2 * z **
                       2 + az4 * z**4)
            y_r = y * z**2 * (by2 * y**2 + bz2 * z**2)
            z_r = y * z**2 * (cy2 * y**2 + cz2 * z**2 +
                              cyz * y*z)
            return self.curve(x_r, y_r, z_r)

        def glv(self, k):
            """GLV scalar multiplication `k`*`self`.

            WARNING: this does not work when k = r for example!!!!!!!
            Reference:
            https://www.iacr.org/archive/crypto2001/21390189.pdf

            """
            k = k % self.curve.r
            if k == 0:
                return self.curve(0, 1, 1)
            m1, m2, m3 = self.curve.glv_params
            m1 = int(m1)
            m2 = int(m2)
            m3 = int(m3)

            b = [floor(mpq(k*m1, self.curve.r)),
                 floor(mpq(k*m2, self.curve.r))]
            s0 = k-b[0] * m1 - b[1] * m3
            s1 = -b[0] * m2 - b[1] * -m1

            p0, p1 = self, self.φ()

            # Multi scalar multiplication with small scalars
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
            prec = [[res, p1], [p0, p0+p1]]
            s0 = int(s0)
            s1 = int(s1)
            n = s0.bit_length()
            while n > 0:
                res = res.dbl()
                n -= 1
                res = res + prec[(s0 >> n) & 1][(s1 >> n) & 1]
                s0 &= ~(1 << n)
                s1 &= ~(1 << n)
            return res

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
