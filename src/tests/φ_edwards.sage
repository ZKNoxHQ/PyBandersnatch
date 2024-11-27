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
mg2ed_x = X/Y
mg2ed_y = (X-1)/(X+1)
ed2mg_x = (1+Y) / (1-Y)
ed2mg_y = (1+Y)/(1-Y)/X
assert ed2mg_x(mg2ed_x, mg2ed_y) == X
assert mg2ed_y(ed2mg_x, ed2mg_y) == Y

# Conversion Weierstrass <-> (twisted) Edwards
ws2ed_x = mg2ed_x(ws2mg_x, ws2mg_y)
ws2ed_y = mg2ed_y(ws2mg_x, ws2mg_y)
ed2ws_x = mg2ws_x(ed2mg_x, ed2mg_y)
ed2ws_y = mg2ws_y(ed2mg_x, ed2mg_y)
assert ws2ed_x(ed2ws_x, ed2ws_y) == X
assert ed2ws_y(ws2ed_x, ws2ed_y) == Y

# Endomorphism in Montgomery model
# (X,Y) -> r_ed_XY, s_ed_XY
r_ed_XY = ws2ed_x(rX(x=ed2ws_x), ed2ws_y * sX(x=ed2ws_x))
s_ed_XY = ws2ed_y(x=rX(x=ed2ws_x))

# Point in Weierstrass model
P = E.random_point()
xP = P[0]
yP = P[1]
assert yP**2 == xP**3 + a*xP + b

# Point in Edwards model
xPP = ws2ed_x(xP, yP)
yPP = ws2ed_y(xP, yP)
assert a_ed * xPP**2 + yPP**2 == 1 + d_ed * xPP**2 * yPP**2

# Endomorphism in Weierstrass model
φxP = rX(x=xP)
φyP = yP * sX(xP, yP)
assert φyP**2 == φxP**3 + a*φxP + b

# Endomorphism in Edwards model
xt = r_ed_XY(xPP, yPP)
yt = s_ed_XY(xPP, yPP)
assert a_ed * xt**2 + yt**2 == 1 + d_ed * xt**2*yt**2

# Consistency
assert ws2ed_x(φxP, φyP) == r_ed_XY(xPP, yPP)
assert ws2ed_y(φxP, φyP) == s_ed_XY(xPP, yPP)


# Efficient endomorphism
# r_ed_XY = α * (x/y)  * (y-y1) * (y-y2)
r_num_y = r_ed_XY.numerator()(x=1).factor()
r_den_y = r_ed_XY.denominator()(x=1).factor()
y1 = r_num_y[0][0].coefficients()[1]
y2 = r_num_y[1][0].coefficients()[1]
α = r_num_y.unit()/r_den_y.unit()
assert α * x/y * (y-y1) * (y-y2) == r_ed_XY

print("α = {}".format(hex(α)))
print("y1 = {}".format(hex(y1)))
print("y2 = {}".format(hex(y2)))


s_num_y = s_ed_XY.numerator()(x=1).factor()
s_den_y = s_ed_XY.denominator()(x=1).factor()
y3 = s_num_y[0][0].coefficients()[1]
y4 = s_num_y[1][0].coefficients()[1]
y5 = s_den_y[0][0].coefficients()[1]
y6 = s_den_y[1][0].coefficients()[1]
β = s_num_y.unit()/s_den_y.unit()
assert β * (y-y3) * (y-y4) / ((y-y5) * (y-y6)) == s_ed_XY

print("β = {}".format(hex(β)))
print("y3 = {}".format(hex(y3)))
print("y4 = {}".format(hex(y4)))
print("y5 = {}".format(hex(y5)))
print("y6 = {}".format(hex(y6)))
