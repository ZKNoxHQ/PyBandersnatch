# -*- coding: utf-8 -*- 
from fp import Fp
from gmpy2 import mpq

class Curve:
    def __init__(self, a, b, r, h):
        # TODO check a and d are defined on the same field?
        self.field = a.field
        self.a = a
        self.b = b
        self.r = r
        self.h = h
        
    def __repr__(self):
        return "Montgomery curve defined by {}*y^2 = x^3 + {}*x^2 + x".format(self.b,self.a)

    __str__ = __repr__

    def __call__(self, x, z):
        return self.Point(x, z, self)
    
    def random(self) :
        x = self.field.random()
        while not(((x**3 + self.a*x**2 + x) / self.b).is_square()):
            x = self.field.random()
        return self.Point(x,1, self)
    
    def j_inv(self):
        return 256 * (self.a**2 - 3)**3 / (self.a**2-4)

    def weierstrass(self) :
        return [-(self.a**2)/3+1, self.a*((self.a**2)*2/9-1)/3]

    # def getPointFromWeierstrass(self, P) :
    #     if P[2] == 0 :
    #         return point.Point(1, 0, self)
    #     if P[2] == 1 :
    #         return point.Point(P[0] - self.a/3, 1, self)
    #     X, Y, Z = P[0]/P[2], P[1]/P[2], 1
    #     return self.getPointFromWeierstrass([X,Y,Z])

    class Point:
        def __init__(self, x, z, curve):
            self.x = x
            self.z = z
            self.curve = curve
        
        def __repr__(self):
            return "Point ({}, {})".format(self.x, self.z)
        
        def __eq__(self,other) :
            return self.x * other.z == self.z * other.x

        def normalize(self) :
            if self.z == 0 :
                return self.curve(1, 0)
            return self.curve(self.x/self.z, 1)

        def in_curve(self, twist=False) :
            if self.z == 0 :
                return True
            x = self.x/self.z
            return ((x**3 + self.curve.a*x**2 + x)/self.curve.b).is_square() + twist == 1

        def weierstrass(self) :
            x = self.x
            z = self.z
            if z == 0 :
                return [0,1,0]
            xn = x/z
            return xn+self.curve.a/3

        def equals(self, other) :
            return (self.z == 0 and other.z == 0) or self.x * other.z == other.x * self.z

        def dbl(self) :
            # [2]P, in the xz-model
            # eprint 2017/212 algo 2
            x = self.x
            z = self.z
            v1 = x+z
            v1 = v1**2
            v2 = x-z
            v2 = v2**2
            xR = v1*v2
            v1 = v1-v2
            v3 = ((self.curve.a+2)/4)*v1
            v3 = v3+v2
            zR = v1*v3
            return self.curve(xR, zR)

        def add(self, Q, PmQ) :
            # P+Q, in the xz-model
            xP, zP = self.x, self.z
            xQ, zQ = Q.x, Q.z
            xm, zm = PmQ.x, PmQ.z
            v0 = xP + zP
            v1 = xQ - zQ
            v1 = v1 * v0
            v0 = xP - zP
            v2 = xQ + zQ
            v2 = v2*v0
            v3 = v1+v2
            v3 = v3**2
            v4 = v1-v2
            v4 = v4**2
            xR = zm *v3
            zR = xm * v4
            return self.curve(xR, zR)

        def naive_mul(self, k) :
            # R = [k]P, in the xz-model
            if k == 0 :
                return self.curve(1, 0) # Smith notation
            k = abs(k) # function does not care about sign(k)
            R0 = self
            R1 = R0.dbl()
            i=0
            R1mR0 = R0
            k_bits = [int(bit) for bit in bin(k)[2:]]
            for i in range(1, len(k_bits)) :
                if k_bits[i] == 0 :
                    [R0, R1] = [R0.dbl(), R0.add(R1, R1mR0)]
                else :
                    [R0, R1] = [R0.add(R1, R1mR0), R1.dbl()]
            if R0.z == 0 :
                return self.curve(1, 0)
            return R0
        
        def multi_scalar_mul(self,k1,other, k2, other_minus_self):
            # Algorithm 9 of https://eprint.iacr.org/2017/212.pdf
            s0, s1, P0, P1, Pm = k1, k2, self, other, other_minus_self

            # TODO can be done faster as we look at points up to a sign, but the condition `if s1<s0` then needs to be considered differently
            while s0<0:
                s0 += self.curve.r
            while s1<0:
                s1 += self.curve.r

            while s0 != 0:
                if s1<s0:
                    s0, s1, P0, P1, Pm, = s1,    s0,         P1,            P0,       Pm
                if s1 <= 4*s0:
                    s0, s1, P0, P1, Pm, = s0,    s1-s0,      P1.add(P0,Pm), P1,       P0
                elif s0 % 2 == s1 % 2:
                    s0, s1, P0, P1, Pm, = s0,    (s1-s0)>>1, P1.add(P0,Pm), P1.dbl(), Pm
                elif s1 % 2 == 0:
                    s0, s1, P0, P1, Pm, = s0,    s1>>1,      P0,            P1.dbl(), P1.add(Pm,P0)
                else:
                    s0, s1, P0, P1, Pm, = s0>>1, s1,         P0.dbl(),      P1,       P0.add(Pm,P1)
            while s1 % 2 == 0:
                s1, P1 = s1>>2, P1.dbl()
            if s1 > 1:
                P1 = s1*P1
            return P1
        
        def is_prime_order_point(self, N) :
            # True if self is of order N, False else.
            return (N*self).z == 0 and self.z != 0

        def φ(self):
            # from eprint 2021/1152
            x = self.x
            z = self.z
            c = self.curve.a+2# TODO can be hardcoded but it is not a big deal
            return self.curve(-(x-z)**2-c*x*z, 2*x*z)
        
        def φ_minus_one(self):
            # see the derivation in `φ.sage`
            α = self.curve.field(13017314467421381532402061398313046228820690393386411611562176812113295071440)
            β = self.curve.field(14989411347484419666605643019079533103863186413725217032868654387860539633484)
            γ = self.curve.field(39953720565912266872856944794434720047230584117801669040511822283402326025498)
            return self.curve(α* self.x * (self.x+ self.z * β)**2, self.z * (self.x + γ * self.z)**2)
        
        def glv(self,k):
            # TODO explain where this comes from
            M1 = [113482231691339203864511368254957623327,   10741319382058138887739339959866629956]
            M2 = [21482638764116277775478679919733259912, -113482231691339203864511368254957623327]
            # inverse*det first row
            N1 = [113482231691339203864511368254957623327, 10741319382058138887739339959866629956]
            det = 13108968793781547619861935127046491459309155893440570251786403306729687672801
            b = [round(mpq(k*N1[0],det)), round(mpq(k*N1[1],det))]
            k1 = k-b[0] * M1[0] - b[1] * M2[0]
            k2 = -b[0] * M1[1] -b[1] * M2[1]
            return self.multi_scalar_mul(k1, self.φ(), k2, self.φ_minus_one())

        def __rmul__(self,k):
            return self.glv(k)