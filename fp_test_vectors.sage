# Bandersnatch, Montgomery model
p = 52435875175126190479447740508185965837690552500527637822603658699938581184513
Fp = GF(p)
a = Fp(0x323f6ab2be9df186b32d48320a504c8841b51474108995013ca91c0b6dfaeb04)
b = Fp(0x49f3b288d4250a92c8a6483fc259aed1cfd09151814f598f415b3f41c93c2850)
a_plus_b = a+b
a_minus_b = a-b
a_mul_b = a*b
a_div_b = a/b
nsq = 1
while Fp(nsq).is_square():
    nsq = -nsq
    if nsq > 0:
        nsq += 1
sqrt_b = b.sqrt()


def test_vector_scalar(k, name):
    print("test_vectors['{}'] = F({})".format(name, hex(k)))


test_vector_scalar(a, 'a')
test_vector_scalar(b, 'b')
test_vector_scalar(a+b, 'a_plus_b')
test_vector_scalar(a*b, 'a_mul_b')
test_vector_scalar(a/b, 'a_div_b')
test_vector_scalar(nsq, 'nsq')
test_vector_scalar(sqrt_b, 'sqrt_b')
