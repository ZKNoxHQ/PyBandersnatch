import unittest
from gmpy2 import random_state, mpz_random, invert, mpz, powmod

class Fp():
    
    def __init__(self, p):
        self.p = mpz(p)
        self._rand_state = random_state()

        # 2-adicity
        self.two_adicity = 0
        while (p-1) % 2**self.two_adicity == 0:
            self.two_adicity += 1
        self.two_adicity -= 1

        # Quadratic non-residue
        self.non_square = mpz(1)
        while self.is_square(self.non_square):
            self.non_square *= -1
            if self.non_square > 0:
                self.non_square += 1

    def __str__(self):
        return "Finite field of characteristic {}".format(p)

    def random(self):
        # probably not secure
        return mpz_random(self._rand_state, self.p)
        
    def is_square(self,x):
        # return 0 if x=0, 1 if x is a non-zero square mod p, and -1 if it is not a square
        # TODO optimize the case (p,q) = (q,p) * something
        if x==0:
            return True
        else:
            return powmod(mpz(x),(self.p-1)//2, self.p) == 1
    
    def mul(self, a, b):
        return (a*b)%self.p
    
    def square(self, a):
        return self.mul(a,a) # TODO can be more efficient ?
    
    def sqrt(self, x):
        # Return the square root (mod p) or None if no solution exists
        # from https://www.lvzl.fr/teaching/2020-21/docs/AA-devoir-1-sujet.pdf
        # and https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm
        if x == 0 :
            return mpz(0)
        if self.p % 4 == 3:
            # TODO ASSERT IS SQUARE BEFORE RETURN
            return powmod(mpz(x),(self.p+1)//4, self.p)
        # Tonelli-Shanks
        z = self.non_square
        M = self.two_adicity
        Q = (self.p-1)//2**M

        c = powmod(z, Q, self.p)
        t = powmod(mpz(x), Q, self.p)
        r = powmod(mpz(x), (Q+1)//2, self.p)
        # TODO REMOVE THIS CHECK
        assert(powmod(t, 1<<self.two_adicity,self.p) == 1)

        while t != 1:
            # Find the least i (0 < i < m) such that t^(2^i) ≡ 1 (mod p)
            t2i = t
            i = 0
            while t2i != 1:
                t2i = self.square(t2i)
                i+=1
            # Update variables for the next iteration
            b = powmod(c, 1<<(M - i - 1), self.p)
            M = i
            c = self.square(b)
            t = self.mul(t,c)
            r = self.mul(r,b)
            
        return r if t == 1 else None
    
    def div(self, x, y):
        x = mpz(x)
        y = mpz(y)
        return self.mul(x,invert(y, self.p))

class TestFp(unittest.TestCase, Fp):

    def test_random(self):
        F = Fp(52435875175126190479447740508185965837690552500527637822603658699938581184513)
        a = F.random()
        b = F.random()
        # a=b with probability 1/p, that is small for a large p
        self.assertFalse(a==b)

    def test_mul(self):
        F = Fp(52435875175126190479447740508185965837690552500527637822603658699938581184513)
        self.assertEqual(6, F.mul(2,3))

    def test_div(self):
        F = Fp(52435875175126190479447740508185965837690552500527637822603658699938581184513)
        x = 1234567890
        y = 9876543210
        expected_division = 35275400988426393170843532131382492504976292327056669892310128774672708665364 # from sage
        self.assertEqual(F.div(x,y), expected_division)

    def test_sqrt(self):
        F = Fp(52435875175126190479447740508185965837690552500527637822603658699938581184513)
        # find a square
        sq = F.random()
        while not(F.is_square(sq)):
            sq = F.random()
        # square-root computation
        root = F.sqrt(sq)
        self.assertEqual(F.square(root), sq)

    def test_is_square(self):
        F = Fp(52435875175126190479447740508185965837690552500527637822603658699938581184513)
        self.assertFalse(F.is_square(F.non_square))
        for i in range(F.non_square):
            self.assertTrue(F.is_square(i))
        
if __name__ == '__main__':
    unittest.main()
