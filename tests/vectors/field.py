# File generated using `sage sage/field.sage > tests/vectors/field.py`.
from src.field import Field
p = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
F = Field(p)
test_vectors = {}
test_vectors['a'] = F(0x323f6ab2be9df186b32d48320a504c8841b51474108995013ca91c0b6dfaeb04)
test_vectors['b'] = F(0x49f3b288d4250a92c8a6483fc259aed1cfd09151814f598f415b3f41c93c2850)
test_vectors['a_plus_b'] = F(0x84575e869257ed14899b869c3082354bdc801c291da92917e045b4e37371353)
test_vectors['a_mul_b'] = F(0x13cbd79f7463c1ca7cfd8c662c6b69d2ef3b2bf383be658c3217117a624973f3)
test_vectors['a_div_b'] = F(0x3ddbe37f9605b5d0f48eabfb2bd0fbd1acd8100d6f7d1172171795192e50fc75)
test_vectors['non_square'] = F(0x5)
test_vectors['sqrt_b'] = F(0x2ca0f83ebf679a313a66801928cad698850f9b0f59a7e7f344e4929eeae32bc7)
