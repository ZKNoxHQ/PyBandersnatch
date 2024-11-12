# GLV in xonly requires P, φ(P) and φ(P)-P
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

Fpx = Fp['x']
x = Fpx.gen()
FpX = Fpx.fraction_field()
X=FpX.gen()

# Endomorphism in Weierstrass model
# (X,Y) -> rX, Y * sX
# Isogeny
alpha = E.division_polynomial(2).roots()[0][0]
P = E.lift_x(alpha)
phi0,phi1 = E.isogeny(P)
E2 = E.isogeny_codomain(P)
# Isomorphism
u = list(set((x**4 - E.a4()/E2.a4()).roots()) & set((x**6 - E.a6()/E2.a6()).roots()))[1][0]
assert u**4 == E.a4()/E2.a4() and u**6 == E.a6()/E2.a6()
rX = phi0(x=X,y=1) * u**2
sX = phi1(x=X,y=1) * u**3
assert sX**2 * (X**3 + a*X + b) == rX**3 + a * rX + b

# Endomorphism in Montgomery model
# (X,Y) -> rmX, Y * smX
# Conversion Montgomery <-> Weierstrass
mg2ws = B*(X+A/3)
ws2mg = X/B-A/3
assert ws2mg(mg2ws) == X
assert mg2ws(ws2mg) == X
# Composition with rX and sX
rmX = ws2mg(rX(mg2ws))
smX = sX(mg2ws)
assert smX**2 * (X**3+A*X**2+X) == rmX**3 +A*rmX**2 + rmX

# Rational map of φ-1 from https://www.iacr.org/archive/eurocrypt2014/84410275/84410275.pdf (page 8)
# The difference x(P-Q) is (xP*xQ * (2*A + xQ + xP) + 2*yP*yQ*B + xP + xQ) /  (xP-xQ)**2
# We set:
# x(P) = X,
# x(Q) = rmX
# y(P) = Y,
# y(Q) = smX * Y
# and use the fact that y(P)*y(Q) = smX * (X**3+A*X**2+X) in order to obtain:
x_φ_minus_one = (X*rmX * (2*A + rmX + X) + 2*smX * (X**3+A*X**2+X) + X + rmX) /  (X-rmX)**2

# The final x(φ-1) = (αt(t+β)², (t+γ)²)
α = Fp(13017314467421381532402061398313046228820690393386411611562176812113295071440)
β = Fp(14989411347484419666605643019079533103863186413725217032868654387860539633484)
γ = Fp(39953720565912266872856944794434720047230584117801669040511822283402326025498)
assert(x_φ_minus_one == α*X*(X+β)**2/(X+γ)**2)
# in projective x-only coordinates: (X,Z) -> (αX(X+Zβ)², Z(X+γZ)²)

print(rmX.numerator().factor())
print(rmX.denominator().factor())
print(smX.numerator().factor())
assert smX.denominator() == rmX.denominator()**2