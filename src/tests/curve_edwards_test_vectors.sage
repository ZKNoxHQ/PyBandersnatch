# Bandersnatch, Montgomery model
Fp = GF(52435875175126190479447740508185965837690552500527637822603658699938581184513)
A = Fp(0x4247698f4e32ad45a293959b4ca17afa4a2d2317e4c6ce5023e1fd63d1b5de98)
B = Fp(5)
r = 13108968793781547619861935127046491459309155893440570251786403306729687672801
h = 4

# Bandersnatch, Weierstrass model
a = B**2 * (1-A**2/3)
b = B**3 * A * (2*A**2/9 - 1) / 3
E = EllipticCurve(Fp, [a, b])
assert (E.j_invariant() == 8000)
assert (E.order() % r == 0)

# Bandersnatch, Edwards model
# 2017-212, equation (6)
a_ed, d_ed = (A+2)/B, (A-2)/B

# generator
gx = 6
while not ((1-a_ed*gx**2)/(1-d_ed*gx**2)).is_square():
    gx = -gx
    if gx > 0:
        gx += 1
gy = ((1-a_ed*gx**2)/(1-d_ed*gx**2)).sqrt()


def to_mg(P):
    u = P[0]/P[2]
    v = P[1]/P[2]
    x = u/B - A/3
    y = v/B**2
    return x, y


def to_ed(P):
    # 2017-212, equations (4) and (5)
    x, y = to_mg(P)
    u = x/y
    v = (x-1)/(x+1)
    return u, v


p = E(31511963179209183026886029814959507395230513391536014203721350106469568871776,
      45347120062487836513813256222005391829075297965413648488198604153937949600247)
q = E(20252373884274151187306374916054971403178400027240398097441693977206289492028,
      51776867226593987565156122032412653740536539297269422625045303422905294395870)
k = 11997154529596648729624281997554038960651754640906483911385998427296165917073
k1 = -45894336995428141233269187797940484884
k2 = 8683555061824981937504960049179714114
λ = -8913659658109529928382530854484400854125314752504019737736543920008458395397

u, v = to_ed(p)
assert a_ed * u**2 + v**2 == 1+d_ed * u**2*v**2


def φ(p):
    # see rmX and smX in `φ.sage`
    x_p, y_p = to_ed(p)
    α = 0x50281ac0f92fc1b29d2a646fe1f5beb21ec0cb08e81f589296d082245cf9382d
    y1 = 0x2123b4c7a71956a2d149cacda650bd7d2516918bf263672811f0feb1e8daef4d
    y2 = 0x52c9f28b828426a561f00d3a63511a882ea712770d9af4d6ee0f014d172510b4
    β = 0x52c9f28b828426a561f00d3a63511a882ea712770d9af4d6ee0f014d172510b4
    y3 = 0x231d3dfa72b4af2d818763ac3629f247733f1d83fd9d3d50e3f6732dae64a8b7
    y4 = 0x50d06958b6e8ce1ab1b2745bd377e5bde07e867f02611eae1c098cd1519b574a
    y5 = 0x150f478323e54389c182cabeeae1fa155d29c5c24b3e595703e518e7e336084c
    y6 = 0x5ede5fd005b839be71b70d491ebfddeff693de40b4c002a7fc1ae7171cc9f7b5
    return α * x_p/y_p * (y_p-y1) * (y_p-y2),  β * (y_p-y3) * (y_p-y4) / ((y_p-y5) * (y_p-y6))

# GENERATION OF TEST VECTORS


def test_vector_point(p, name, ws=True):
    if ws:
        [x, y] = to_ed(p)
    else:
        [x, y] = p
    print("test_vectors['{}'] = E(F({}), F({}), F(1))".format(
        name, hex(x), hex(y)))


def test_vector_scalar(k, name):
    print("test_vectors['{}'] = {}".format(name, hex(k)))


print("F = Field({})".format(hex(Fp.characteristic())))
print("a = F({})".format(a_ed))
print("d = F({})".format(d_ed))
print("r = 0x1cfb69d4ca675f520cce760202687600ff8f87007419047174fd06b52876e7e1")
print("h = 4")
print("E = CurveEdwards(a, d, r, h)")
print("test_vectors = {}")
test_vector_point(p, 'p')
test_vector_point(q, 'q')
test_vector_point(2*p, 'p_double')
test_vector_point(p+q, 'p_plus_q')
test_vector_point(p-q, 'p_minus_q')
test_vector_point(λ*p, 'φ_p')
test_vector_point(k*p, 'k_times_p')
test_vector_point(k1*p, 'k1_times_p')
test_vector_point(k2*p, 'k2_times_p')
test_vector_point(k1*p + k2*q, 'k1_times_p_plus_k2_times_q')
test_vector_scalar(k, 'k')
test_vector_scalar(k1, 'k1')
test_vector_scalar(k2, 'k2')
test_vector_scalar(λ, 'λ')
