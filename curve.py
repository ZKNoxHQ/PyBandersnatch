# -*- coding: utf-8 -*- 
import unittest
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

        def __rmul__(self, k) :
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
            s0, s1, P0, P1, Pm = abs(k1), abs(k2), self, other, other_minus_self

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
            x = self.x
            z = self.z
            c = self.curve.a+2#field(0x4247698f4e32ad45a293959b4ca17afa4a2d2317e4c6ce5023e1fd63d1b5de9a) 
            return self.curve(-(x-z)**2-c*x*z, 2*x*z)
        
        def φ_minus_one(self):
            α = self.curve.field(13017314467421381532402061398313046228820690393386411611562176812113295071440)
            β = self.curve.field(14989411347484419666605643019079533103863186413725217032868654387860539633484)
            γ = self.curve.field(39953720565912266872856944794434720047230584117801669040511822283402326025498)
            return self.curve(α* self.x * (self.x+ self.z * β)**2, self.z * (self.x + γ * self.z)**2)
        
        def glv(self,k):
            # N0 and det obtained using:
            # ```
            # M = Matrix([[-L,1], [r,0]]).LLL()
            # det = det(M)
            # N = det * M.inv()
            # N0 = [N[0][0], N[0][1]]
            # ```
            M1 = [113482231691339203864511368254957623327,   10741319382058138887739339959866629956]
            M2 = [21482638764116277775478679919733259912, -113482231691339203864511368254957623327]
            # inverse*det first row
            N1 = [113482231691339203864511368254957623327, 10741319382058138887739339959866629956]
            det = 13108968793781547619861935127046491459309155893440570251786403306729687672801
            b = [round(mpq(k*N1[0],det)), round(mpq(k*N1[1],det))]
            k1 = k-b[0] * M1[0] - b[1] * M2[0]
            k2 = -b[0] * M1[1] -b[1] * M2[1]
            print("k ={}".format(k))
            print("k1={}".format(k1))
            print("k2={}".format(k2))
            return self.multi_scalar_mul(k1, self.φ(), k2, self.φ_minus_one())

class TestCurve(unittest.TestCase):

    def set_up_curve(self):
        # bandersnatch curve in montgomery model
        F = Fp(52435875175126190479447740508185965837690552500527637822603658699938581184513)
        a = F(0x4247698f4e32ad45a293959b4ca17afa4a2d2317e4c6ce5023e1fd63d1b5de98)
        b = F(0x300c3385d13bedb7c9e229e185c4ce8b1dd3b71366bb97c30855c0aa41d62727)
        r = 13108968793781547619861935127046491459309155893440570251786403306729687672801 
        h = 4
        return Curve(a, b, r, h)
    
    def test_j_invariant(self):
        E = self.set_up_curve()
        self.assertEqual(E.j_inv(), 8000)

    def test_random(self):
        E = self.set_up_curve()
        for i in range(10):
            P = E.h*E.random()
            self.assertTrue(P.is_prime_order_point(E.r))

    def test_φ_norm(self):
        E = self.set_up_curve()
        P = E.h*E.random()
        φ2P = P.φ().φ()
        self.assertEqual(φ2P, -2*P)

    def test_φ_eigenvalue(self):
        E = self.set_up_curve()
        F = E.field
        xP = F(12082691309842786344904703102623539864205261229594144938888324074174613840300)
        P = E(xP,1)
        self.assertTrue(P.is_prime_order_point(E.r))
        λ = 8913659658109529928382530854484400854125314752504019737736543920008458395397
        self.assertEqual(P.φ(), λ*P)

    def test_φ(self):
        E = self.set_up_curve()
        F = E.field
        xP = F(12082691309842786344904703102623539864205261229594144938888324074174613840300)
        P = E(xP,1)
        φP = P.φ()
        self.assertEqual(φP, E(F(11458922273495097112533829568725487588121633355211186266238294705768134023326),1))

    def test_φ_minus_one(self):
        E = self.set_up_curve()
        F = E.field
        xP = F(33458486013926512409896242145198845310109221311939503997839007609527653194379)
        P = E(xP, 1)
        φ_minus_one_P = P.φ_minus_one()
        expected_x = F(24387954571941077226850081312365069165536914211052705257835068138877697985688)
        self.assertEqual(φ_minus_one_P, E(expected_x, 1))

    def test_φ_minus_one_2(self):
        E = self.set_up_curve()
        F = E.field
        xP = F(12082691309842786344904703102623539864205261229594144938888324074174613840300)
        P = E(xP, 1)
        φ_minus_one_P = P.φ_minus_one()
        expected_x = F(24033430638574135559207661009508680256329241072539898145614847227501177452871)
        self.assertEqual(φ_minus_one_P, E(expected_x, 1))

    def test_add(self):
        E = self.set_up_curve()
        F = E.field
        xP = F(33458486013926512409896242145198845310109221311939503997839007609527653194379)
        xQ = F(28482916637127167529646451618185459303123407338943306674359938301245425470092)
        xPmQ = F(32798664934949466701208494782160952539180648511133632801776750414975815504392)
        xPpQ = F(12304079527202110574548882275124160605549627239042994825351407972975127435146)

        P = E(xP,1)
        Q = E(xQ,1)
        PmQ = E(xPmQ, 1)
        PpQ = P.add(Q,PmQ)
        self.assertEqual(PpQ, E(xPpQ,1))

    def test_dbl(self):
        E = self.set_up_curve()
        F = E.field
        xP = F(33458486013926512409896242145198845310109221311939503997839007609527653194379)
        x2P = F(12082691309842786344904703102623539864205261229594144938888324074174613840300)

        P = E(xP,1)
        Pdbl = P.dbl()
        self.assertEqual(Pdbl, E(x2P,1))

    def test_scalar_mul(self):
        E = self.set_up_curve()
        F = E.field
        xP = F(33458486013926512409896242145198845310109221311939503997839007609527653194379)
        k = 2345676543234567654
        xkP = F(38783438532518312076050495050645672027147987262410212490358095175313247827184)

        P = E(xP,1)
        kP = k*P
        self.assertEqual(kP, E(xkP,1))

        xQ = F(28482916637127167529646451618185459303123407338943306674359938301245425470092)
        l = 234567890987659876987658765
        xlQ = F(39150136521723522337504363078544378164314913244247808907542546539324382480286)

        Q = E(xQ,1)
        lQ = l*Q
        self.assertEqual(lQ, E(xlQ, 1))

    def test_multi_scalar_mul(self):
        E = self.set_up_curve()
        F = E.field
        xP = F(33458486013926512409896242145198845310109221311939503997839007609527653194379)
        xQ = F(28482916637127167529646451618185459303123407338943306674359938301245425470092)
        xPmQ = F(32798664934949466701208494782160952539180648511133632801776750414975815504392)
        xk1Pk2Q = F(39864008317233774338731736742297798582452053800324681545427967023194223256838)

        P = E(xP,1)
        Q = E(xQ,1)
        PmQ = E(xPmQ, 1)
        k1 = 2345676543234567654
        k2 = 234567890987659876987658765
        k1_p_plus_k2_q = P.multi_scalar_mul(k1,Q,k2, PmQ)
        self.assertEqual(k1_p_plus_k2_q, E(xk1Pk2Q,1))

    def test_multi_scalar_mul_2(self):
        E = self.set_up_curve()
        F = E.field
        xP = F(12082691309842786344904703102623539864205261229594144938888324074174613840300)
        P = E(xP,1)
        Q = P.φ()
        xPmQ = F(24033430638574135559207661009508680256329241072539898145614847227501177452871)
        PmQ = E(xPmQ,1)

        k1 = -45894336995428141233269187797940484884
        k2 = 8683555061824981937504960049179714114
        k1_p_plus_k2_q = P.multi_scalar_mul(k1,Q,k2, PmQ)
        self.assertEqual(k1_p_plus_k2_q, E(F(24025662955270513540262306395182467712802082813859915548603039872586808886210),1))

    # def test_glv(self):
    #     E = self.set_up_curve()
    #     F = E.field        
    #     xP = F(12082691309842786344904703102623539864205261229594144938888324074174613840300)
    #     xkP =F(24025662955270513540262306395182467712802082813859915548603039872586808886210)
    #     k = 11997154529596648729624281997554038960651754640906483911385998427296165917073

    #     P = E(xP,1)
    #     kP = -45894336995428141233269187797940484884*P + 8683555061824981937504960049179714114 * P.φ()
    #     self.assertEqual(kP, E(xkP,1))
    
    # def test_bench(self):
    #     from time import time
    #     # from bandersnatch paper (montgomery)
    #     F = Fp(52435875175126190479447740508185965837690552500527637822603658699938581184513)
    #     a = F(0x4247698f4e32ad45a293959b4ca17afa4a2d2317e4c6ce5023e1fd63d1b5de98)
    #     b = F(0x300c3385d13bedb7c9e229e185c4ce8b1dd3b71366bb97c30855c0aa41d62727)
    #     r = 13108968793781547619861935127046491459309155893440570251786403306729687672801 
    #     h = 4
    #     E = Curve(a,b,r,h)
        
    #     xP = F(12082691309842786344904703102623539864205261229594144938888324074174613840300)
    #     P = E(xP,1)
    #     k = 1310896879378154761986193512704649000000155893440570251786403306729687672801

    #     t = time()
    #     for i in range(100):
    #         kP_glv = P.glv(k)
    #     print(time()-t)

    #     t = time()
    #     for i in range(100):
    #         kP_naive = k*P
    #     print(time()-t)
    #     self.assertEqual(kP_glv, kP_naive)
    

        
if __name__ == '__main__':
    unittest.main()
