# File generated using `sage sage/montgomery.sage > tests/vectors/montgomery.py`.
from src.field import Field
from src.curve.montgomery25519 import Montgomery25519
F = Field(57896044618658097711785492504343953926634992332820282019728792003956564819949)
a = F(486662)
b = F(1)
r = 0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3ed
h = 8
E = Montgomery25519(a, b, r, h)
test_vectors = {}
test_vectors['p'] = E(F(0x460dfa666eef17ce354700e223e2ee569b63ec0f5891533013fb47e56296e239), 1)
test_vectors['q'] = E(F(0x7800aac5e16d2479479a7739cb68d5b15702f1326051ef5ea1be2f105f652fe8), 1)
test_vectors['p_double'] = E(F(0x4d1b82d9a06919aa419d08d0ec50d52f8810f1dfab41b1ece85c615c8a341ca3), 1)
test_vectors['p_plus_q'] = E(F(0x362e1a8d8652803a93ae839b19cbbc20fe8678fa6152beb5f61cea53dda899dc), 1)
test_vectors['p_minus_q'] = E(F(0x6e92ddfdc0da00180cefb68a4943ad61d3dd32d0ebd033b83572df7b7c693a35), 1)
test_vectors['k_times_p'] = E(F(0x442f1439e6efb73033bccda19f040b3d25067c543364ff9da01da62900ff657a), 1)
test_vectors['k1_times_p'] = E(F(0x1820dddf208738f9c30f25b21262fcf70a1cd090f0b30718ecc2ffb6f326c158), 1)
test_vectors['k2_times_p'] = E(F(0x6883197fa524088b9043ea28ace64767ffd63113c9da15f43b5827da4eceee2d), 1)
test_vectors['k1_times_p_plus_k2_times_q'] = E(F(0xe3899d8a71d12a554a2384d305aa387722218cd3bcf28908839cc6e451f9501), 1)
test_vectors['k'] = 0x1a862619b8224e61eb24bb583c84ce04913064d37308623924c7a64fcdc9f191
test_vectors['k1'] = -0x2286ed83a0b1545d1b7788921e40bb14
test_vectors['k2'] = 0x6886451b4aa55294c626bb34d42e242
