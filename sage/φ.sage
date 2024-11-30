# Montgomery curve coefficients of Bandersnatch
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


# Bandersnatch, Edwards model
# 2017-212, equation (6)
a_ed, d_ed = (A+2)/B, (A-2)/B
# Isomorphism to get a=5
assert (a_ed/-5).is_square()
sqrt_a_ed_5 = (a_ed/-5).sqrt()
a_ed, d_ed = Fp(-5), d_ed / a_ed * -5


Fpxy = Fp['x', 'y']
x, y = Fpxy.gens()
FpXY = Fpxy.fraction_field()
X, Y = FpXY.gens()

# Endomorphism in Weierstrass model
# (X,Y) -> rX, Y * sX
# Isogeny
alpha = E.division_polynomial(2).roots()[0][0]
P = E.lift_x(alpha)
phi0, phi1 = E.isogeny(P)
E2 = E.isogeny_codomain(P)
# Isomorphism
u = (E.a4()/E2.a4()).sqrt().sqrt()
assert u**4 == E.a4()/E2.a4() and u**6 == E.a6()/E2.a6()
rX = phi0(X, 1) * u**2
sX = phi1(X, 1) * u**3
assert sX**2 * (X**3 + a*X + b) == rX**3 + a * rX + b

# Conversion Montgomery <-> Weierstrass
mg2ws_x = B*(X+A/3)
mg2ws_y = B**2 * Y
ws2mg_x = X/B-A/3
ws2mg_y = Y/B**2
assert ws2mg_x(mg2ws_x, mg2ws_y) == X
assert mg2ws_y(ws2mg_x, ws2mg_y) == Y
# Conversion Montgomery <-> (twisted) Edwards
mg2ed_x = X/Y * sqrt_a_ed_5
mg2ed_y = (X-1)/(X+1)
ed2mg_x = (1+Y) / (1-Y)
ed2mg_y = (1+Y)/(1-Y) / (X/sqrt_a_ed_5)
assert ed2mg_x(mg2ed_x, mg2ed_y) == X
assert mg2ed_y(ed2mg_x, ed2mg_y) == Y
# Conversion Weierstrass <-> (twisted) Edwards
ws2ed_x = mg2ed_x(ws2mg_x, ws2mg_y)
ws2ed_y = mg2ed_y(ws2mg_x, ws2mg_y)
ed2ws_x = mg2ws_x(ed2mg_x, ed2mg_y)
ed2ws_y = mg2ws_y(ed2mg_x, ed2mg_y)
assert ws2ed_x(ed2ws_x, ed2ws_y) == X
assert ed2ws_y(ws2ed_x, ws2ed_y) == Y

# Endomorphism in Edwards model
# (X,Y) -> r_ed_XY, s_ed_XY
r_ed_XY = ws2ed_x(rX(x=ed2ws_x), ed2ws_y * sX(x=ed2ws_x))
s_ed_XY = ws2ed_y(x=rX(x=ed2ws_x))

# Endomorphism in Montgomery model
r_mg_XY = ws2mg_x(rX(x=mg2ws_x), mg2ws_y * sX(x=mg2ws_x))
s_mg_XY = ws2mg_y(y=mg2ws_y * sX(x=mg2ws_x))


# P and φ(P) in Weierstrass model
P = E(31511963179209183026886029814959507395230513391536014203721350106469568871776,
      45347120062487836513813256222005391829075297965413648488198604153937949600247)
x_ws_P = P[0]
y_ws_P = P[1]
assert y_ws_P**2 == x_ws_P**3 + a*x_ws_P + b
x_ws_φ_P = rX(x=x_ws_P)
y_ws_φ_P = y_ws_P * sX(x_ws_P, y_ws_P)
assert y_ws_φ_P**2 == x_ws_φ_P**3 + a*x_ws_φ_P + b

# P and φ(P) in Montgomery model
x_mg_P = ws2mg_x(x_ws_P, y_ws_P)
y_mg_P = ws2mg_y(x_ws_P, y_ws_P)
assert B*y_mg_P**2 == x_mg_P**3 + A * x_mg_P**2 + x_mg_P
x_mg_φ_P = r_mg_XY(x=x_mg_P, y=y_mg_P)
y_mg_φ_P = s_mg_XY(x_mg_P, y_mg_P)
assert B * y_mg_φ_P**2 == x_mg_φ_P**3 + A*x_mg_φ_P**2 + x_mg_φ_P

# P and φ(P) in Edwards model
x_ed_P = mg2ed_x(x_mg_P, y_mg_P)
y_ed_P = mg2ed_y(x_mg_P, y_mg_P)
assert a_ed * x_ed_P**2 + y_ed_P**2 == 1 + d_ed * x_ed_P**2 * y_ed_P**2
x_ed_φ_P = r_ed_XY(x_ed_P, y_ed_P)
y_ed_φ_P = s_ed_XY(x_ed_P, y_ed_P)
assert a_ed * x_ed_φ_P**2 + y_ed_φ_P**2 == 1 + d_ed * x_ed_φ_P**2*y_ed_φ_P**2

# Consistency
assert ws2ed_x(x_ws_φ_P, y_ws_φ_P) == r_ed_XY(x_ed_P, y_ed_P)
assert ws2ed_y(x_ws_φ_P, y_ws_φ_P) == s_ed_XY(x_ed_P, y_ed_P)
assert ws2mg_x(x_ws_φ_P, y_ws_φ_P) == r_mg_XY(x_mg_P, y_mg_P)
assert ws2mg_y(x_ws_φ_P, y_ws_φ_P) == s_mg_XY(x_mg_P, y_mg_P)


# Projective version of the endomorphism in Edwards model
Fpxyz = Fp['x', 'y', 'z']
FpXYZ = FractionField(Fpxyz)
X, Y, Z = FpXYZ.gens()
# (x:y:z) -> (x/z,y/z,1) -> (r_ed_XY(x/z), y/z s_ed_XY(x/z), 1)
rXYZ = r_ed_XY(x=X/Z, y=Y/Z)
sXYZ = s_ed_XY(x=X/Z, y=Y/Z)
# (x:y:z) -> (z * a(x,z) : y * b(x,z) : z * c(x,z))
a = rXYZ.numerator() * sXYZ.denominator()
b = sXYZ.numerator() * rXYZ.denominator()
c = rXYZ.denominator() * sXYZ.denominator()

# Hand-made simplification
a_yz = a(x=1)
b_yz = Fpxyz(b(x=1) / (Y*Z**2))
c_yz = Fpxyz(c(x=1) / (Y*Z**2))
assert a == X * a_yz
assert b == Y * Z**2 * b_yz
assert c == Y * Z**2 * c_yz

[ay4, ay2z2, az4] = a_yz.coefficients()
[by2, bz2] = b_yz.coefficients()
[cy2, cz2] = c_yz.coefficients()
print("# Edwards model endomorphism")
print("ay4 = {}".format(hex(ay4)))
print("ay2z2 = {}".format(hex(ay2z2)))
print("az4 = {}".format(hex(az4)))
print("by2 = {}".format(hex(by2)))
print("bz2 = {}".format(hex(bz2)))
print("cy2 = {}".format(hex(cy2)))
print("cz2 = {}".format(hex(cz2)))
print()

assert (x_ed_P * a_yz(y=y_ed_P, z=1) /
        (y_ed_P * 1**2 * c_yz(y=y_ed_P, z=1)) == x_ed_φ_P)
assert (y_ed_P * 1**2 * b_yz(y=y_ed_P, z=1) /
        (y_ed_P * 1**2 * c_yz(y=y_ed_P, z=1)) == y_ed_φ_P)

'''
# Rational map of φ-1 from https://www.iacr.org/archive/eurocrypt2014/84410275/84410275.pdf (page 8)
# The difference x(P-Q) is (xP*xQ * (2*A + xQ + xP) + 2*yP*yQ*B + xP + xQ) /  (xP-xQ)**2
# We set:
# x(P) = X, x(Q) = rmX, y(P) = Y, y(Q) = smX * Y
# and use the fact that y(P)*y(Q) = smX * (X**3+A*X**2+X) in order to obtain:
x_φ_minus_one = (X*r_mg_XY * (2*A + r_mg_XY + X) + 2*s_mg_XY *
                 (X**3+A*X**2+X) + X + r_mg_XY) / (X-r_mg_XY)**2
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
'''

# GLV decomposition
λ = -8913659658109529928382530854484400854125314752504019737736543920008458395397
[M1, M2] = Matrix([[-λ, 1], [r, 0]]).LLL()
assert ((M1[0] + λ * M1[1]) % r == 0)
N1 = (r * Matrix([M1, M2])**-1)[0]
print("# GLV parameters")
print("M1 = [{}, {}]".format(M1[0], M1[1]))
print("M2 = [{}, {}]".format(M2[0], M2[1]))
print("N1 = [{}, {}]".format(N1[0], N1[1]))
