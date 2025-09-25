# Bandersnatch, Montgomery model
p = (1 << 255)-19
Fp = GF(p)
a = Fp(0xd467b32c86a5f4340f65400e372747f0f5a5b1493266e1e173578e9bacb9b23)
b = Fp(0x7d9e36316caca29cf9525c784ba072eebecb2f7607b47f26e8ca588f8bea0b3c)
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


print("# File generated using `sage sage/ed25519_field.sage > tests/vectors/ed25519_field.py`.")
print("from src.field import Field")
print("p = {}".format(hex(p)))
print("F = Field(p)")
print("test_vectors = {}")
test_vector_scalar(a, 'a')
test_vector_scalar(b, 'b')
test_vector_scalar(a+b, 'a_plus_b')
test_vector_scalar(a*b, 'a_mul_b')
test_vector_scalar(a/b, 'a_div_b')
test_vector_scalar(nsq, 'non_square')
test_vector_scalar(sqrt_b, 'sqrt_b')
