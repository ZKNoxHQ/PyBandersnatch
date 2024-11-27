# -*- coding: utf-8 -*-
import argparse
from field import Field
from gmpy2 import mpq, mpz, mod


def constant_time_swap(swap_flag, a, b):
    """Constant time swap function."""
    swap_flag = int(bool(swap_flag))
    # TODO check if this mask is okay
    mask = - swap_flag

    a_new = (a & ~mask) | (b & mask)
    b_new = (b & ~mask) | (a & mask)
    return a_new, b_new


class CurveEdwards:
    def __init__(self, a, d, r, h):
        self.field = a.field
        self.a = a
        self.d = d
        self.r = r
        self.h = h
        self.a24 = (self.a+2)/4
        self.generator = self.Point(self.field(1), self.field(
            0x31574b8e6115fd25314d9139bc32a29822b510278f5c8a72d6a3d41aa253bc6), self.field(1), self)

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

            Double-and-add algorithm.
            TODO not constant time

            """
            if k == 0:
                return self.curve(0, 1, 0)
            if k < 0:
                k = -k
                self = self.neg()
            res = self.curve(0, 1, 0)
            temp = self
            k_bits = [int(bit) for bit in bin(k)[-1:1:-1]]
            for i in range(len(k_bits)):
                if k_bits[i] == 1:
                    res = res.add(temp)
                temp = temp.dbl()
            if res.z == 0:
                return self.curve(1, 0, 0)
            return res

        def multi_scalar_mul(self, k1, other, k2):
            """Multi scalar multiplication `k1` * `self` + `k2` * `other`.

            TODO not constant time.

            """
            s0, s1, p0, p1 = k1, k2, self, other

            if s0 < 0:
                s0 = -s0
                p0 = p0.neg()
            if s1 < 0:
                s1 = -s1
                p1 = p1.neg()

            if s0 == 0 and s1 == 0:
                return self.curve(0, 1, 0)

            if s1 > s0:
                s0, p0, s1, p1 = s1, p1, s0, p0
            # s1 ≤ s0

            res = self.curve(0, 1, 0)
            temp = [[res, p0], [p1, p0.add(p1)]]
            s0 = int(s0)
            s1 = int(s1)
            n = s0.bit_length()
            print(s1.bit_length())
            while s0 > 0:
                res = res.dbl()
                n -= 1
                print("n={}".format(n))
                print(bin(s0))
                print(bin(s1))
                print("s0>>n = {}".format(s0 >> n))
                print("s1>>n = {}".format(s1 >> n))
                res = res.add(temp[s0 >> n][s1 >> n])
                s0 ^= (1 << n)
                s1 ^= (1 << n)
            return res

       # def is_prime_order(self, N):
       #     """Returns the boolean corresponding to `self.order() == N`."""
       #     return self.naive_mul(N).z == 0 and self.z != 0

       # def φ(self):
       #     """Endomorphism sqrt(-2).

       #     Reference:
       #     https://eprint.iacr.org/2021/1152.pdf page 6.

       #     """
       #     x = self.x
       #     z = self.z
       #     c = self.curve.a+2  # TODO can be hard-coded
       #     return self.curve(-(x-z)**2-c*x*z, 2*x*z)

       # def φ_minus_one(self):
       #     """Endomorphism sqrt(-2) - [1].

       #     More information in the file `φ.sage`.

       #     """
       #     α = self.curve.field(
       #         13017314467421381532402061398313046228820690393386411611562176812113295071440)
       #     β = self.curve.field(
       #         14989411347484419666605643019079533103863186413725217032868654387860539633484)
       #     γ = self.curve.field(
       #         39953720565912266872856944794434720047230584117801669040511822283402326025498)
       #     return self.curve(α * self.x * (self.x + self.z * β)**2, self.z * (self.x + γ * self.z)**2)

       # def glv(self, k, constant_time=False):
       #     """GLV scalar multiplication `k`*`self`.

       #     A constant time option is available.
       #     Reference:
       #     https://www.iacr.org/archive/crypto2001/21390189.pdf
       #     More information in the file `φ.sage`.

       #     """
       #     if k == 0:
       #         return self.curve(1, 0)
       #     M1 = [113482231691339203864511368254957623327,
       #           10741319382058138887739339959866629956]
       #     M2 = [21482638764116277775478679919733259912, -
       #           113482231691339203864511368254957623327]
       #     N1 = [113482231691339203864511368254957623327,
       #           10741319382058138887739339959866629956]
       #     b = [round(mpq(k*N1[0], self.curve.r)),
       #          round(mpq(k*N1[1], self.curve.r))]
       #     k1 = k-b[0] * M1[0] - b[1] * M2[0]
       #     k2 = -b[0] * M1[1] - b[1] * M2[1]
       #     return self.multi_scalar_mul(k1, self.φ(), k2, self.φ_minus_one(), constant_time=constant_time)

       # def __rmul__(self, k, constant_time=False):
       #     """Scalar multiplication with the scalar give first.

       #     Computed using GLV.

       #     """
       #     return self.glv(k, constant_time=constant_time)

       # # def slow_add(self, q):
       # #     """Compute the addition `self` ± `q`.

       # #     It is not a differential addition, but requires the computation of the y-coordinates with sqrt.
       # #     TODO optimized.
       # #     Reference:
       # #     https://www.iacr.org/archive/eurocrypt2014/84410275/84410275.pdf page 8.

       # #     """
       # #     x_p = self.x/self.z
       # #     x_q = q.x/q.z
       # #     a = self.curve.a
       # #     b = self.curve.b

       # #     y_p = ((x_p**3 + a*x_p**2 + x_p)/b).sqrt()
       # #     y_q = ((x_q**3 + a*x_q**2 + x_q)/b).sqrt()

       # #     return self.curve(b * (x_q * y_p - x_p*y_q)**2 / (x_p*x_q*(x_p-x_q)**2), 1)
