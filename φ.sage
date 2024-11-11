# GLV in xonly requires P, φ(P) and φ(P)-P
# Computed from https://www.iacr.org/archive/eurocrypt2014/84410275/84410275.pdf (page 8)

p = 52435875175126190479447740508185965837690552500527637822603658699938581184513
Fp = GF(p)
A = Fp(0x4247698f4e32ad45a293959b4ca17afa4a2d2317e4c6ce5023e1fd63d1b5de98)
B = Fp(0x300c3385d13bedb7c9e229e185c4ce8b1dd3b71366bb97c30855c0aa41d62727)
r = 13108968793781547619861935127046491459309155893440570251786403306729687672801 
h = 4

# Conversion to Weierstrass
a = B**2 * (1-A**2/3)
b = B**3 * A * (2*A**2/9 - 1) / 3 
E = EllipticCurve(Fp, [a,b])
assert(E.j_invariant() == 8000)
assert(E.order()%r == 0)

# 1. φ in eprint 2021/1152 (weierstrass model y² = x³+ax+b)
aa = Fp(-3763200000)
bb = Fp(-78675968000000)

Fpt = Fp['t']
t = Fpt.gen()
u = (t**4 - aa/a).roots()[0][0]

t0 = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfefffffffef10be001
u0 = 0x50281ac0f92fc1b20fd897a16bf2b9e132bdcb06721c589296cf82245cf9382d
Ft = Fpt.fraction_field()
T=Ft.gen()

r1t = u0**2 * (T**2  + 44800 * T + 2257920000)
r2t = T + 44800
r3t = u0**3 * (T**2  + 2*44800 * T + t0)
E2 = EllipticCurve(Fp, [aa,bb])
P = 4 * E2.random_point()
assert P.order() == r
φP = E2(r1t(P[0])/r2t(P[0]), P[1] * r3t(P[0])/r2t(P[0])**2)

# In the first weierstrass model that is equivalent to our montgomery model
rT, sT = r1t(u**2*T)/u**2 / r2t(u**2*T), r3t(u**2*T)/r2t(u**2*T)**2
Q = 4 * E.random_point()
assert Q.order()==r
φQ = E(rT(Q[0]), Q[1] * sT(Q[0]))
φ2Q = E(rT(φQ[0]), φQ[1] * sT(φQ[0]))
assert φ2Q == -2*Q
assert sT**2 * (T**3 + a*T + b) == rT**3 + a * rT + b
# (x,y) -> r(x), y s(x) endomorphism of E

# In Montgomery model
# MG2WS (x,y) -> (B(x+A/3), B²y)
# WS2MG (u,v) -> (u/B-A/3, v/B²)
mg2ws = B*(T+A/3)
ws2mg = T/B-A/3
assert ws2mg(mg2ws) == T
assert mg2ws(ws2mg) == T

r2T = ws2mg(rT(mg2ws))
s2T = sT(mg2ws)
assert s2T**2 * (T**3+A*T**2+T) == r2T**3 +A*r2T**2 + r2T
# (x,y) -> r2T(x), y * s2T(x) endomorphism of E in montgomery model

# Rational map of φ-1
# The difference x(P-Q) is (xP*xQ * (2*A + xQ + xP) + 2*yP*yQ*B + xP + xQ) /  (xP-xQ)**2
# Using x(P) = T and x(Q) = r2T(T)  and y(Q) = s2T * y defined above, we obtain:
x_φ_minus_one = (T*r2T * (2*A + r2T + T) + 2*s2T * (T**3+A*T**2+T) + T + r2T) /  (T-r2T)**2

# The final x(φ-1) = (αt(t+β)², (t+γ)²)
α = Fp(13017314467421381532402061398313046228820690393386411611562176812113295071440)
β = Fp(14989411347484419666605643019079533103863186413725217032868654387860539633484)
γ = Fp(39953720565912266872856944794434720047230584117801669040511822283402326025498)
assert(x_φ_minus_one == α*T*(T+β)**2/(T+γ)**2)
# in projective x-only coordinates: (x,z) mapsto (αx(x+zβ)², z(x+γz)²)
