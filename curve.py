# -*- coding: utf-8 -*-
import argparse
from field import Field
from gmpy2 import mpq, mpz


def constant_time_swap(swap_flag, a, b):
    """Constant time swap function."""
    swap_flag = int(bool(swap_flag))
    # TODO check if this mask is okay
    mask = - swap_flag

    a_new = (a & ~mask) | (b & mask)
    b_new = (b & ~mask) | (a & mask)
    return a_new, b_new


class Curve:
    def __init__(self, a, b, r, h):
        self.field = a.field
        self.a = a
        self.b = b
        self.r = r
        self.h = h
        self.a24 = (self.a+2)/4
        self.generator = self.Point(self.field(0xa), 1, self)

    def __repr__(self):
        return "Montgomery curve defined by {}*y^2 = x^3 + {}*x^2 + x".format(self.b, self.a)

    __str__ = __repr__

    def __call__(self, x, z):
        if isinstance(x, int) or isinstance(x, mpz):
            x = self.field(x)
        if isinstance(z, int) or isinstance(z, mpz):
            z = self.field(z)
        return self.Point(x, z, self)

    def random(self):
        """Returns a random point of `self`."""
        x = self.field.random()
        while not (((x**3 + self.a*x**2 + x) / self.b).is_square()):
            x = self.field.random()
        return self.Point(x, self.field(1), self)

    def j_inv(self):
        """Returns the j-invariant of `self`.

        Reference:
        https://eprint.iacr.org/2017/212.pdf page 3.

        """
        return 256 * (self.a**2 - 3)**3 / (self.a**2-4)

    def weierstrass(self):
        """Weierstrass coefficients of `self`.

        Reference:
        https://eprint.iacr.org/2017/212.pdf page 6.

        """
        return [-(self.a**2)/3+1, self.a*((self.a**2)*2/9-1)/3]

    class Point:
        def __init__(self, x, z, curve):
            self.x = x
            self.z = z
            self.curve = curve

        def __repr__(self):
            return "Point ({}, {})".format(self.x, self.z)

        def __eq__(self, other):
            """Return the equality boolean between `self` and `other`.

            Points are represented in projective coordinates.
            The equality is defined modulo {±1}.
            """
            return self.x * other.z == other.x * self.z

        def normalize(self):
            """Affine representation of the projective point."""
            if self.z == 0:
                return self.curve(1, 0)
            return self.curve(self.x/self.z, 1)

        def in_curve(self, twist=False):
            """Returns the curve membership boolean."""
            if self.z == 0:
                return True
            x = self.x/self.z
            # `not()` because self.curve.b is a non-square!
            return not ((x**3 + self.curve.a*x**2 + x).is_square()) + twist == 1

        def dbl(self):
            """Doubling algorithm.

            Reference:
            https://eprint.iacr.org/2017/212.pdf algorithm 2.

            """
            x = self.x
            z = self.z
            v1 = x+z
            v1 = v1**2
            v2 = x-z
            v2 = v2**2
            x_r = v1*v2
            v1 = v1-v2
            v3 = self.curve.a24*v1
            v3 = v3+v2
            z_r = v1*v3
            return self.curve(x_r, z_r)

        def add(self, q, p_minus_q):
            """Differential addition algorithm.

            Requires p_minus_q.
            Reference:
            https://eprint.iacr.org/2017/212.pdf algorithm 1.

            """
            x_p, z_p = self.x, self.z
            x_q, zQ = q.x, q.z
            xm, zm = p_minus_q.x, p_minus_q.z
            v0 = x_p + z_p
            v1 = x_q - zQ
            v1 = v1 * v0
            v0 = x_p - z_p
            v2 = x_q + zQ
            v2 = v2*v0
            v3 = v1+v2
            v3 = v3**2
            v4 = v1-v2
            v4 = v4**2
            x_r = zm * v3
            z_r = xm * v4
            return self.curve(x_r, z_r)

        def naive_mul(self, k):
            """Scalar multiplication `k` * `self`.

            Reference:
            https://eprint.iacr.org/2017/212.pdf algorithm 3.

            """
            if k == 0:
                return self.curve(1, 0)
            k = abs(k)  # computation modulo {±1}
            r0 = self
            r1 = r0.dbl()
            i = 0
            r1_minus_r0 = r0
            k_bits = [int(bit) for bit in bin(k)[2:]]
            for i in range(1, len(k_bits)):
                if k_bits[i] == 0:
                    [r0, r1] = [r0.dbl(), r0.add(r1, r1_minus_r0)]
                else:
                    [r0, r1] = [r0.add(r1, r1_minus_r0), r1.dbl()]
            if r0.z == 0:
                return self.curve(1, 0)
            return r0

        def constant_time_point_swap(self, other, swap_flag):
            p_1, p_2 = self, other
            p_1_x, p_2_x = constant_time_swap(
                swap_flag, p_1.x.value, p_2.x.value)
            p_1_z, p_2_z = constant_time_swap(
                swap_flag, p_1.z.value, p_2.z.value)
            return self.curve(p_1_x, p_1_z), self.curve(p_2_x, p_2_z)

        def mul_rfc_7748(self, k):
            """Scalar multiplication `k` * `self` following RFC 7748.

            Reference:
            https://datatracker.ietf.org/doc/html/rfc7748
            """
            if k == 0:
                return self.curve(1, 0)
            k = abs(k)  # computation modulo {±1}
            r0 = self
            r1 = r0.dbl()
            i = 0
            r1_minus_r0 = r0
            k_bits = [int(bit) for bit in bin(k)[2:]]
            swap = 0
            for i in range(1, len(k_bits)):
                swap ^= k_bits[i]
                r0, r1 = r0.constant_time_point_swap(r1, swap)
                r0, r1 = r0.dbl(), r0.add(r1, r1_minus_r0)
            # TODO is it `not(swap)`` or `swap`? Papers say `swap`` but here it works with `not(swap)`
            r0, r1 = r0.constant_time_point_swap(r1, not (swap))
            if r0.z == 0:
                return self.curve(1, 0)
            return r0

        def multi_scalar_mul(self, k1, other, k2, other_minus_self, constant_time=False):
            """Multi scalar multiplication `k1` * `self` + `k2` * `other`.

            Reference:
            https://eprint.iacr.org/2017/212.pdf Algorithm 9.

            """
            s0, s1, p0, p1, pm = k1, k2, self, other, other_minus_self

            # TODO is it constant-time ? is swapping constant time really important?
            if s0 < 0:
                s0 = -s0
                pm = p0.add(p1, pm)
            if s1 < 0:
                s1 = -s1
                pm = p0.add(p1, pm)

            if s0 == 0 and s1 == 0:
                return self.curve(1, 0)

            while s0 != 0:
                if s1 < s0:
                    if constant_time:
                        _ = p0.dbl()
                        _ = p0.add(p1, pm)
                    s0, s1, p0, p1, pm, = s1, s0,  p1,  p0, pm
                if s1 <= 4*s0:
                    if constant_time:
                        _ = p0.dbl()
                    s0, s1, p0, p1, pm, = s0, s1 - s0,  p1.add(p0, pm), p1, p0
                elif s0 % 2 == s1 % 2:
                    s0, s1, p0, p1, pm, = s0, (s1 -
                                               s0) >> 1, p0.add(p1, pm), p1.dbl(), pm
                elif s1 % 2 == 0:
                    s0, s1, p0, p1, pm, = s0, s1 >> 1,  p0,  p1.dbl(), p1.add(pm, p0)
                else:
                    s0, s1, p0, p1, pm, = s0 >> 1, s1,  p0.dbl(), p1, p0.add(pm, p1)
            while s1 % 2 == 0:
                if constant_time:
                    _ = p0.dbl()
                s1, p1 = s1 >> 1, p1.dbl()
            if s1 > 1:
                p1 = p1.naive_mul(s1)
            return p1

        def is_prime_order(self, N):
            """Returns the boolean corresponding to `self.order() == N`."""
            return self.naive_mul(N).z == 0 and self.z != 0

        def φ(self):
            """Endomorphism sqrt(-2).

            Reference:
            https://eprint.iacr.org/2021/1152.pdf page 6.

            """
            x = self.x
            z = self.z
            c = self.curve.a+2  # TODO can be hard-coded
            return self.curve(-(x-z)**2-c*x*z, 2*x*z)

        def φ_minus_one(self):
            """Endomorphism sqrt(-2) - [1].

            More information in the file `φ.sage`.

            """
            α = self.curve.field(
                13017314467421381532402061398313046228820690393386411611562176812113295071440)
            β = self.curve.field(
                14989411347484419666605643019079533103863186413725217032868654387860539633484)
            γ = self.curve.field(
                39953720565912266872856944794434720047230584117801669040511822283402326025498)
            return self.curve(α * self.x * (self.x + self.z * β)**2, self.z * (self.x + γ * self.z)**2)

        def glv(self, k, constant_time=False):
            """GLV scalar multiplication `k`*`self`.

            A constant time option is available.
            Reference:
            https://www.iacr.org/archive/crypto2001/21390189.pdf
            More information in the file `φ.sage`.

            """
            if k == 0:
                return self.curve(1, 0)
            M1 = [113482231691339203864511368254957623327,
                  10741319382058138887739339959866629956]
            M2 = [21482638764116277775478679919733259912, -
                  113482231691339203864511368254957623327]
            N1 = [113482231691339203864511368254957623327,
                  10741319382058138887739339959866629956]
            b = [round(mpq(k*N1[0], self.curve.r)),
                 round(mpq(k*N1[1], self.curve.r))]
            k1 = k-b[0] * M1[0] - b[1] * M2[0]
            k2 = -b[0] * M1[1] - b[1] * M2[1]
            return self.multi_scalar_mul(k1, self.φ(), k2, self.φ_minus_one(), constant_time=constant_time)

        def __rmul__(self, k, constant_time=False):
            """Scalar multiplication with the scalar give first.

            Computed using GLV.

            """
            return self.glv(k, constant_time=constant_time)

        # def slow_add(self, q):
        #     """Compute the addition `self` ± `q`.

        #     It is not a differential addition, but requires the computation of the y-coordinates with sqrt.
        #     TODO optimized.
        #     Reference:
        #     https://www.iacr.org/archive/eurocrypt2014/84410275/84410275.pdf page 8.

        #     """
        #     x_p = self.x/self.z
        #     x_q = q.x/q.z
        #     a = self.curve.a
        #     b = self.curve.b

        #     y_p = ((x_p**3 + a*x_p**2 + x_p)/b).sqrt()
        #     y_q = ((x_q**3 + a*x_q**2 + x_q)/b).sqrt()

        #     return self.curve(b * (x_q * y_p - x_p*y_q)**2 / (x_p*x_q*(x_p-x_q)**2), 1)
