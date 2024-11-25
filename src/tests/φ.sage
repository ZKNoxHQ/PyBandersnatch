# GLV in xonly requires P, φ(P) and φ(P)-P
p = 52435875175126190479447740508185965837690552500527637822603658699938581184513
Fp = GF(p)
A = Fp(0x4247698f4e32ad45a293959b4ca17afa4a2d2317e4c6ce5023e1fd63d1b5de98)
B = Fp(5)
r = 13108968793781547619861935127046491459309155893440570251786403306729687672801
h = 4

# Conversion to Weierstrass
a = B**2 * (1-A**2/3)
b = B**3 * A * (2*A**2/9 - 1) / 3
E = EllipticCurve(Fp, [a, b])
assert (E.j_invariant() == 8000)
assert (E.order() % r == 0)

Fpx = Fp['x']
x = Fpx.gen()
FpX = Fpx.fraction_field()
X = FpX.gen()

# Endomorphism in Weierstrass model
# (X,Y) -> rX, Y * sX
# Isogeny
alpha = E.division_polynomial(2).roots()[0][0]
P = E.lift_x(alpha)
phi0, phi1 = E.isogeny(P)
E2 = E.isogeny_codomain(P)
# Isomorphism
u = list(set((x**4 - E.a4()/E2.a4()).roots()) &
         set((x**6 - E.a6()/E2.a6()).roots()))[1][0]
assert u**4 == E.a4()/E2.a4() and u**6 == E.a6()/E2.a6()
rX = phi0(x=X, y=1) * u**2
sX = phi1(x=X, y=1) * u**3
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
assert smX**2 * (X**3+A*X**2+X) == rmX**3 + A*rmX**2 + rmX

# Rational map of φ-1 from https://www.iacr.org/archive/eurocrypt2014/84410275/84410275.pdf (page 8)
# The difference x(P-Q) is (xP*xQ * (2*A + xQ + xP) + 2*yP*yQ*B + xP + xQ) /  (xP-xQ)**2
# We set:
# x(P) = X, x(Q) = rmX, y(P) = Y, y(Q) = smX * Y
# and use the fact that y(P)*y(Q) = smX * (X**3+A*X**2+X) in order to obtain:
x_φ_minus_one = (X*rmX * (2*A + rmX + X) + 2*smX *
                 (X**3+A*X**2+X) + X + rmX) / (X-rmX)**2

(x_φ_minus_one_numerator_1, mul_1), (x_φ_minus_one_numerator_2,
                                     mul_2) = x_φ_minus_one.numerator().factor()
(x_φ_minus_one_denominator, mul_3), = x_φ_minus_one.denominator().factor()

# The final x(φ-1) = (αt(t+β)², (t+γ)²)
α = x_φ_minus_one.numerator().factor().unit()
β = x_φ_minus_one_numerator_2.list()[0]
γ = x_φ_minus_one_denominator.list()[0]
assert (x_φ_minus_one == α*X*(X+β)**2/(X+γ)**2)

# in projective x-only coordinates: (X,Z) -> (αX(X+Zβ)², Z(X+γZ)²)
Fpxz = Fp['x', 'z']
FpXZ = FractionField(Fpxz)
X, Z = FpXZ.gens()
x_φ_minus_one_numerator_1_xz = FpXZ(x_φ_minus_one_numerator_1)(x=X/Z)
x_φ_minus_one_numerator_2_xz = FpXZ(x_φ_minus_one_numerator_2)(x=X/Z)
x_φ_minus_one_denominator_xz = FpXZ(x_φ_minus_one_denominator)(x=X/Z)

# GLV decomposition
λ = 8913659658109529928382530854484400854125314752504019737736543920008458395397
[M1, M2] = Matrix([[-λ, 1], [r, 0]]).LLL()
assert ((M1[0] + λ * M1[1]) % r == 0)
N1 = (r * Matrix([M1, M2])**-1)[0]
