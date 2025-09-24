# File generated using `sage sage/ed25519_field.sage > tests/vectors/ed25519_field.py`.
from src.field import Field
p = 0x7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffed
F = Field(p)
test_vectors = {}
test_vectors['a'] = F(
    0xd467b32c86a5f4340f65400e372747f0f5a5b1493266e1e173578e9bacb9b23)
test_vectors['b'] = F(
    0x7d9e36316caca29cf9525c784ba072eebecb2f7607b47f26e8ca588f8bea0b3c)
test_vectors['a_plus_b'] = F(
    0xae4b164351701e03a48b0792f12e76dce258a8a9adaed44ffffd17946b5a672)
test_vectors['a_mul_b'] = F(
    0xba39586a8adcc76d0ffdf873def2a9fc9c44f096059cae903d273569e3f58fd)
test_vectors['a_div_b'] = F(
    0x426a24066c2229fe334b81046e277cce0f05a058be37ca724f488fb44f6ff89f)
test_vectors['non_square'] = F(0x2)
test_vectors['sqrt_b'] = F(
    0x25d8af83cba27469fea9d53ad9bea15e92c0634e578cef56a2ca39fb4dae4771)
