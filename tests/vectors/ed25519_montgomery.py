# File generated from RFC 7748.
from src.curve.edwards import Edwards
from src.field import Field
from src.curve.montgomery import Montgomery
F = Field(2**255-19)
a = F(486662)
b = F(1)
r = 7237005577332262213973186563042994240857116359379907606001950938285454250989
h = 8
E = Montgomery(a, b, r, h, glv=False)
test_vectors = {}
test_vectors['p'] = E(
    F(0xe6db6867583030db3594c1a424b15f7c726624ec26b3353b10a903a6d0ab1c4c), 1)
