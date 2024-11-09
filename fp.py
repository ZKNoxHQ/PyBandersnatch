import unittest
from gmpy2 import random_state, mpz_random, invert, mpz, powmod, sign, mod

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
        self.non_square = self.Element(mpz(1), self)
        while self.non_square.is_square():
            self.non_square += self.Element(1, self)

    def __str__(self):
        return "Finite field of characteristic {}".format(self.p)

    def __call__(self, value):
        return self.Element(value, self)
    
    def random(self):
        # probably not secure
        return self.Element(mpz_random(self._rand_state, self.p), self)
                


    class Element:
        def __init__(self, value, field):
            self.value = value % field.p
            self.field = field
        
        def __eq__(self,other):
            if isinstance(other, self.field.Element) and self.field == other.field:
                return self.value == other.value 
            raise TypeError("Cannot add elements from different fields.")
        
        def __add__(self, other):
            if isinstance(other, self.field.Element) and self.field == other.field:
                return self.field(self.value + other.value) # TODO SOMETIMES DO `% self.field.p`
            raise TypeError("Cannot add elements from different fields.")

        def __neg__(self):
            return self.field(-self.value % self.field.p)
        
        def __sub__(self, other):
            if isinstance(other, self.field.Element) and self.field == other.field:
                return self.field(self.value - other.value) # TODO SOMETIMES DO `% self.field.p)`
            raise TypeError("Cannot subtract elements from different fields.")

        def __mul__(self, other):
            if isinstance(other, self.field.Element) and self.field == other.field:
                return self.field((self.value * other.value) % self.field.p)
            raise TypeError("Cannot multiply elements from different fields.")

        def __pow__(self, exponent):
            # Raise the element to the power of `exponent` modulo `modulus`
            result = powmod(self.value, exponent, self.field.p)
            return self.field(result)

        def is_square(self):
            # return 0 if x=0, 1 if x is a non-zero square mod p, and -1 if it is not a square
            # TODO optimize the case (p,q) = (q,p) * something
            if self.value==0:
                return True
            else:
                return self ** ((self.field.p-1)>>1) == self.field(1)

        def __truediv__(self, other):
            if isinstance(other, self.field.Element) and self.field == other.field:
                inverse_other = invert(other.value, self.field.p)
                return self.field((self.value * inverse_other) % self.field.p)
            # TODO division by zero error?
            raise TypeError("Cannot divide elements from different fields.")

        def __repr__(self):
            return f"{self.value}"
    
        def sqrt(self):
            # Return the square root (mod p) or None if no solution exists
            # from https://www.lvzl.fr/teaching/2020-21/docs/AA-devoir-1-sujet.pdf
            # and https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm

            if self.value == 0 :
                return self
            
            if self.field.p % 4 == 3:
                # TODO ASSERT IS SQUARE BEFORE RETURN
                return powmod(self.value,(self.field.p+1)//4, self.field.p)
            
            # Tonelli-Shanks
            z = self.field.non_square
            M = self.field.two_adicity
            Q = (self.field.p-1)//2**M

            c = z**Q
            t = self**Q
            r = self**((Q+1)>>1)
            # TODO REMOVE THIS CHECK
            assert(t**(1<<self.field.two_adicity) == self.field(1))

            while t != self.field(1):
                # Find the least i (0 < i < m) such that t^(2^i) ≡ 1 (mod p)
                t2i = t
                i = 0
                while t2i != self.field(1):
                    t2i = t2i**2
                    i+=1
                    
                # Update variables for the next iteration
                b = c**(1<<(M-i-1))
                M = i
                c = b*b
                t = t*c
                r = r*b
                
            return r if t == self.field(1) else None
    

class TestFp(unittest.TestCase, Fp):

    def test_random(self):
        F = Fp(52435875175126190479447740508185965837690552500527637822603658699938581184513)
        a = F.random()
        b = F.random()
        # a=b with probability 1/p, that is small for a large p
        self.assertFalse(a==b)

    def test_mul(self):
        F = Fp(52435875175126190479447740508185965837690552500527637822603658699938581184513)
        a = F(12345)
        b = F(54321)
        a_mul_b = F(670592745)
        self.assertEqual(a_mul_b, a*b)

    def test_div(self):
        F = Fp(52435875175126190479447740508185965837690552500527637822603658699938581184513)
        a = F(1234567890)
        b = F(9876543210)
        a_div_b = F(35275400988426393170843532131382492504976292327056669892310128774672708665364) # from sage
        self.assertEqual(a/b, a_div_b)

    def test_sqrt(self):
        F = Fp(52435875175126190479447740508185965837690552500527637822603658699938581184513)
        # find a square
        sq = F.random()
        while not(sq.is_square()):
            sq = F.random()
        # square-root computation
        root = sq.sqrt()#F.sqrt(sq)
        self.assertEqual(root*root, sq)

    def test_is_square(self):
        F = Fp(52435875175126190479447740508185965837690552500527637822603658699938581184513)
        self.assertFalse(F.non_square.is_square())
        for i in range(F.non_square.value):
            self.assertTrue(F(i).is_square())
        
if __name__ == '__main__':
    unittest.main()
