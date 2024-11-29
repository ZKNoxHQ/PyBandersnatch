# File generated using `sage sage/montgomery.sage > tests/vectors/montgomery.py`.
from src.field import Field
from src.curve.montgomery import Montgomery
F = Field(0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001)
a = F(0x4247698f4e32ad45a293959b4ca17afa4a2d2317e4c6ce5023e1fd63d1b5de98)
b = F(5)
r = 0x1cfb69d4ca675f520cce760202687600ff8f87007419047174fd06b52876e7e1
h = 4
E = Montgomery(a, b, r, h)
test_vectors = {}
test_vectors['p'] = E(F(0x3d65d0dd26c89b6c1e111a4b7d4875a8fc60cc1a2191beba4f6c11f2fb1c2cbf), 1)
test_vectors['q'] = E(F(0x4f9ace66d99bf333dd1eb464c902fcc9217372d8902583ff2e65fb3df88fdbb8), 1)
test_vectors['p_double'] = E(F(0x5f6aff9f5f3a3dedf24000ef9899766ccdadc2183675dac057d3d754faf6174), 1)
test_vectors['p_plus_q'] = E(F(0x6770424666db6f35101a26a4ea455b8886346935cb00298962781fa1f5f654af), 1)
test_vectors['p_minus_q'] = E(F(0x45d99562094ef0bfe61ffb151ce57cf76c7c3b247179cb3ebcee99c3e0a59889), 1)
test_vectors['φ_p'] = E(F(0x48850dbbc447be3c1f566275495c29ad981ed87de6818a73e1918a130a9b524d), 1)
test_vectors['φ_minus_one_p'] = E(F(0x1521dba053303cfa98d4e0f358cea947c69d96f5bf8a8d3a075694b4c9b6558d), 1)
test_vectors['k_times_p'] = E(F(0x35d13b10631a46276be4c289a63570de6c552376605b2881e445031c90a96a7e), 1)
test_vectors['k1_times_p'] = E(F(0x367b03140f02455d6d00a7924a65635f717f6354c23ccc2d5ed2b11e4b736c9d), 1)
test_vectors['k2_times_p'] = E(F(0xd6643d2acb959c8f4306b8398ddf8102dbd9fcc52c17817e1bf7943cb8a40c7), 1)
test_vectors['k1_times_p_plus_k2_times_q'] = E(F(0x182a5d0dce9c87ecd9050d911058c4ec30cbc680f94e5c347da06eb068c2dd15), 1)
test_vectors['k'] = 0x1a862619b8224e61eb24bb583c84ce04913064d37308623924c7a64fcdc9f191
test_vectors['k1'] = -0x2286ed83a0b1545d1b7788921e40bb14
test_vectors['k2'] = 0x6886451b4aa55294c626bb34d42e242
test_vectors['λ'] = 0x13b4f3dc4a39a493edf849562b38c72bcfc49db970a5056ed13d21408783df05
