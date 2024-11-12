# -*- coding: utf-8 -*- 
import unittest
from fp import Fp
from curve import Curve

class TestCurve(unittest.TestCase):

    def set_up_curve(self):
        # Bandersnatch curve in montgomery model
        F = Fp(52435875175126190479447740508185965837690552500527637822603658699938581184513)
        a = F(0x4247698f4e32ad45a293959b4ca17afa4a2d2317e4c6ce5023e1fd63d1b5de98)
        b = F(0x300c3385d13bedb7c9e229e185c4ce8b1dd3b71366bb97c30855c0aa41d62727)
        r = 13108968793781547619861935127046491459309155893440570251786403306729687672801 
        h = 4
        E = Curve(a,b,r,h)
        # obtained from `sage curve_test_vectors.sage`
        P = E(F(0x1ab68f993bb722ee7080fb5c9183f65311fd186f6090ec8e135e007723dc61ac), 1)
        Q = E(F(0x3ef8c26fccab7da37bb8d1ef48eca32b9fd8f69770e2dec7e43ea7aa52604a8c), 1)
        Pdbl = E(F(0x458e35834e64c6275fe8d0750ee93e91abf07a4030e437e15cb738beb03a64c3), 1)
        PpQ = E(F(0x1d0b0e9f15b24e5153bae79aa9da7a4332c0cae106f80d673f3c55af4519c4d5), 1)
        PmQ = E(F(0x1703b9d77e39fea52eb7b91b8afd4f032aad2b7e8c7cbb4422882a73ed44365), 1)
        φP = E(F(0x19558529c55e5bbd8e4c0ff655d57c38b40af2eaf62cba0e95bac2ef921b189e), 1)
        φ_minus_one_P = E(F(0x3522708f7799920e9bcd244a44612883230a7e0e01040f0f8858188a3817b947), 1)
        kP = E(F(0x351e0b1814900f6256103c91ce3cf3f668785ffdbd2b5f9a9c94089a0f1103c2), 1)
        k1P = E(F(0xd597b82d17d409eebd4c4658dcbb34bb58c6c9c84aea1fafa6bf65cb98b6082), 1)
        k2P = E(F(0x20d5933448366dd6da2b84d1a060fc7ec2e3e84b6ec44d7a625d93b61b50f40d), 1)
        k1Ppk2Q = E(F(0x40a029d09426fb5d96337e243a42bed0d2eda4c2e1051a587a83d2bac7c2c2fa), 1) 
        k = 0x1a862619b8224e61eb24bb583c84ce04913064d37308623924c7a64fcdc9f191
        k1 = -0x2286ed83a0b1545d1b7788921e40bb14
        k2 = 0x6886451b4aa55294c626bb34d42e242
        λ = 0x13b4f3dc4a39a493edf849562b38c72bcfc49db970a5056ed13d21408783df05
        return E,P, [Q,Pdbl,PpQ,PmQ,φP, φ_minus_one_P, kP, k1P, k2P, k1Ppk2Q, k, k1, k2, λ]
    
    def test_j_invariant(self):
        E,P,test_vecs = self.set_up_curve()
        self.assertEqual(E.j_inv(), 8000)

    def test_random(self):
        E,_,_ = self.set_up_curve()
        for i in range(10):
            P = E.random().naive_mul(E.h)
            self.assertTrue(P.is_prime_order_point(E.r))

    def test_is_prime_order_point(self):
        E,P,test_vecs = self.set_up_curve()
        self.assertTrue(P.is_prime_order_point(E.r))

    def test_φ_norm(self):
        E,P,test_vecs = self.set_up_curve()
        φ2P = P.φ().φ()
        self.assertEqual(φ2P, test_vecs[1])

    def test_φ_eigenvalue(self):
        E,P,test_vecs = self.set_up_curve()
        self.assertEqual(P.φ(), test_vecs[-1]*P)

    def test_φ(self):
        E,P,test_vecs = self.set_up_curve()
        φP = P.φ()
        self.assertEqual(φP, test_vecs[4])

    def test_φ_minus_one(self):
        E,P,test_vecs = self.set_up_curve()
        φ_minus_one_P = P.φ_minus_one()
        self.assertEqual(φ_minus_one_P, test_vecs[5])

    def test_add(self):
        E,P,test_vecs = self.set_up_curve()
        Q = test_vecs[0]
        PmQ = test_vecs[3]
        PpQ = P.add(Q,PmQ)
        self.assertEqual(PpQ, test_vecs[2])

    def test_dbl(self):
        E,P,test_vecs = self.set_up_curve()
        Pdbl = P.dbl()
        self.assertEqual(Pdbl, test_vecs[1])

    def test_scalar_mul(self):
        E,P,test_vecs = self.set_up_curve()
        k = test_vecs[10]
        kP = P.naive_mul(k)
        self.assertEqual(kP, test_vecs[6])

        k1 = test_vecs[11]
        k1P = P.naive_mul(k1)
        self.assertEqual(k1P, test_vecs[7])

        k2 = test_vecs[12]
        k2P = P.naive_mul(k2)
        self.assertEqual(k2P, test_vecs[8])

    def test_multi_scalar_mul(self):
        E,P,test_vecs = self.set_up_curve()
        F = E.field
        Q = test_vecs[0]
        PmQ = test_vecs[3]
        k1 = test_vecs[11]
        k2 = test_vecs[12]
        k1_p_plus_k2_q = P.multi_scalar_mul(k1,Q,k2, PmQ)
        self.assertEqual(k1_p_plus_k2_q, test_vecs[9])

    def test_glv(self):
        E,P,test_vecs = self.set_up_curve()
        k = test_vecs[10]
        kP = P.glv(k)
        self.assertEqual(kP, test_vecs[6])
    
    def test_bench(self):
        from time import time
        E,P,test_vecs = self.set_up_curve()
        k = 1199715452959664211345788999754554038123456678896406483911385998427296165917073
        # test_vecs[10]

        t = time()
        for i in range(200):
            kP_glv = k*P # glv
        time_glv = time() - t

        t = time()
        for i in range(200):
            kP_naive = P.naive_mul(k)
        time_naive = time()-t
        print("GLV is {:.0f}% faster than a scalar multiplication.".format(time_glv/time_naive*100))
        self.assertEqual(kP_glv, kP_naive)
    
if __name__ == '__main__':
    unittest.main()
