Fp = GF((1 << 255)-19)
Fpxyz = Fp['x', 'y', 'z']
FpXYZ = FractionField(Fpxyz)
X, Y, Z = FpXYZ.gens()
r = 13108968793781547619861935127046491459309155893440570251786403306729687672801
h = 4

# Ed25519, Twisted-Edwards model
a_ed, d_ed = Fp(-1), Fp(-121665/121666)

# MODELS OF ELLIPTIC CURVES

# Conversion Montgomery <-> (twisted) Edwards
ed2mg_x = (1+Y) / (1-Y)
ed2mg_y = (1+Y)/(1-Y) / X
mg2ed_x = X/Y
mg2ed_y = (X-1)/(X+1)

a_mg, b_mg = 2*(a_ed+d_ed)/(a_ed-d_ed), 4/(a_ed-d_ed)

x_te, y_te = 15112221349535400772501151409588531511454012693041857206046113283949847762202, 46316835694926478169428394003475163141307993866256225615783033603165251855960
x_mg, y_mg = ed2mg_x(Y=y_te), ed2mg_y(X=x_te, Y=y_te)

assert b_mg * y_mg**2 == x_mg**3 + a_mg * x_mg**2 + x_mg


# # Conversion Montgomery <-> Weierstrass
# mg2ws_x = B*(X+A/3)
# mg2ws_y = B**2 * Y
# ws2mg_x = X/B-A/3
# ws2mg_y = Y/B**2

# # Conversion Weierstrass <-> (twisted) Edwards
# ws2ed_x = mg2ed_x(ws2mg_x, ws2mg_y, Z)
# ws2ed_y = mg2ed_y(ws2mg_x, ws2mg_y, Z)
# ed2ws_x = mg2ws_x(ed2mg_x, ed2mg_y, Z)
# ed2ws_y = mg2ws_y(ed2mg_x, ed2mg_y, Z)

# # generator
# gx = 2
# while not ((1-a_ed*gx**2)/(1-d_ed*gx**2)).is_square():
#     gx = -gx
#     if gx > 0:
#         gx += 1
# gy = ((1-a_ed*gx**2)/(1-d_ed*gx**2)).sqrt()


# def to_mg(P):
#     u = P[0]/P[2]
#     v = P[1]/P[2]
#     x = u/B - A/3
#     y = v/B**2
#     return x, y


# def to_ed(P):
#     # 2017-212, equations (4) and (5)
#     x, y = to_mg(P)
#     u = x/y * sqrt_a_ed_5
#     v = (x-1)/(x+1)
#     return u, v


# k = 11997154529596648729624281997554038960651754640906483911385998427296165917073

# u, v = to_ed(p)
# assert a_ed * u**2 + v**2 == 1+d_ed * u**2*v**2


# def φ(p):
#     # endomorphism in affine coordinates as in 2021/1152.pdf
#     x_p, y_p = p[0]/p[2], p[1]/p[2]
#     alpha = E.division_polynomial(2).roots()[0][0]
#     P = E.lift_x(alpha)
#     phi0, phi1 = E.isogeny(P)
#     E2 = E.isogeny_codomain(P)
#     # Isomorphism
#     u = (E.a4()/E2.a4()).sqrt().sqrt()
#     assert u**4 == E.a4()/E2.a4() and u**6 == E.a6()/E2.a6()
#     rX = phi0(y=1) * u**2
#     sX = phi1(y=1) * u**3
#     Fp = E.base_field()
#     return (Fp(rX(x=x_p, y=y_p)), y_p * Fp(sX(x=x_p, y=y_p)), 1)


# # GENERATION OF TEST VECTORS
# xx, yy, zz = φ(p)
# assert zz == 1 and yy**2 == xx**3 + a*xx + b


# def test_vector_point(p, name, ws=True, projective=False):
#     if ws:
#         [x, y] = to_ed(p)
#         print("test_vectors['{}'] = E(F({}), F({}), F(1))".format(
#             name, hex(x), hex(y)))
#     else:
#         if projective == False:
#             [x, y] = p
#             print("test_vectors['{}'] = E(F({}), F({}), F(1))".format(
#                 name, hex(x), hex(y)))
#         else:
#             [x, y, z] = p
#             print("test_vectors['{}'] = E(F({}), F({}), F({}))".format(
#                 name, hex(x), hex(y), hex(z)))


# def test_vector_scalar(k, name):
#     print("test_vectors['{}'] = {}".format(name, hex(k)))


# print("# File generated using `sage sage/bandersnatch_edwards.sage > tests/vectors/bandersnatch_edwards.py`.")
# print("from src.field import Field")
# print("from src.curve.edwards import Edwards")
# print("F = Field({})".format(hex(Fp.characteristic())))
# print("a = F({})".format(a_ed))
# print("d = F({})".format(d_ed))
# print("r = 0x1cfb69d4ca675f520cce760202687600ff8f87007419047174fd06b52876e7e1")
# print("h = 4")
# print("E = Edwards(a, d, r, h, glv=True)")
# print("test_vectors = {}")
# test_vector_point(p, 'p')
# test_vector_point(q, 'q')
# test_vector_point(2*p, 'p_double')
# test_vector_point(p+q, 'p_plus_q')
# test_vector_point(p-q, 'p_minus_q')
# test_vector_point(φ(p), 'φ_p')
# test_vector_point(k*p, 'k_times_p')
# test_vector_point(k1*p, 'k1_times_p')
# test_vector_point(k2*p, 'k2_times_p')
# test_vector_point(k1*p + k2*q, 'k1_times_p_plus_k2_times_q')
# test_vector_scalar(k, 'k')
# test_vector_scalar(k1, 'k1')
# test_vector_scalar(k2, 'k2')
# test_vector_scalar(λ, 'λ')

# # small order point
# test_vector_point((0, Fp(-1)), "p_order_2_1", False)
# test_vector_point((0, 1, 0), "p_order_2_2", False, True)
# test_vector_point((1, 0, 0), "p_order_2_3", False, True)


# # small x point
# xx = Fp(1)
# yy = sqrt((1-a_ed*xx**2)/(1-d_ed*xx**2))
# xx_ws = 29473314690983384253693538305232603171909485083805813020990978841153575434625
# yy_ws = 44104446522130670290417281820901020888519367242619929491530721381946199062276
# small_p_ws = E(xx_ws, yy_ws)
# assert small_p_ws.order() != r
# test_vector_point((xx, yy), "small_p", False)
# test_vector_point(2 * small_p_ws, "small_p_dbl")
