# Bandersnatch, Montgomery model
Fp = GF((1 << 255)-19)
A = Fp(486662)
B = Fp(1)
r = (1 << 252)+27742317777372353535851937790883648493
h = 8

# Bandersnatch, Weierstrass model
a = B**2 * (1-A**2/3)
b = B**3 * A * (2*A**2/9 - 1) / 3
E = EllipticCurve(Fp, [a, b])
# assert (E.j_invariant() == 8000)
assert (E.order() % r == 0)

# Bandersnatch, Edwards model
# 2017-212, equation (6)
a_ed, d_ed = (A+2)/B, (A-2)/B
# Isomorphism to get a=5
assert (a_ed/-5).is_square()
sqrt_a_ed_5 = (a_ed/-5).sqrt()
a_ed, d_ed = Fp(-5), d_ed/a_ed*-5

# generator
gx = 2
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
    u = x/y * sqrt_a_ed_5
    v = (x-1)/(x+1)
    return u, v


# p = 8*E.random_point()
# print("p = E({}, {})".format(p[0], p[1]))
# q = 8*E.random_point()
# print("q = E({}, {})".format(q[0], q[1]))
p = E(50985278150407432344796787702948612699476491671352921418144363141586907039370,
      15253002477052688986620419946290387236033221949375887382230120676776923013400)
q = E(15681357382615180190462123999934171085226024319595864033686213637078275740748,
      47407999444187020738561942407576495663998743082126352560608723174318508905182)
k = 11997154529596648729624281997554038960651754640906483911385998427296165917073
k1 = -45894336995428141233269187797940484884
k2 = 8683555061824981937504960049179714114

u, v = to_ed(p)
assert a_ed * u**2 + v**2 == 1+d_ed * u**2*v**2


def test_vector_point(p, name, ws=True, projective=False):
    if ws:
        [x, y] = to_ed(p)
        print("test_vectors['{}'] = E(F({}), F({}), F(1))".format(
            name, hex(x), hex(y)))
    else:
        if projective == False:
            [x, y] = p
            print("test_vectors['{}'] = E(F({}), F({}), F(1))".format(
                name, hex(x), hex(y)))
        else:
            [x, y, z] = p
            print("test_vectors['{}'] = E(F({}), F({}), F({}))".format(
                name, hex(x), hex(y), hex(z)))


def test_vector_scalar(k, name):
    print("test_vectors['{}'] = {}".format(name, hex(k)))


print("# File generated using `sage sage/edwards.sage > tests/vectors/edwards.py`.")
print("from src.field import Field")
print("from src.curve.edwards25519 import Edwards25519")
print("F = Field({})".format(hex(Fp.characteristic())))
print("a = F({})".format(a_ed))
print("d = F({})".format(d_ed))
print("r = {}".format(r))
print("h = {}".format(h))
print("E = Edwards25519(a, d, r, h)")
print("test_vectors = {}")
test_vector_point(p, 'p')
test_vector_point(q, 'q')
test_vector_point(2*p, 'p_double')
test_vector_point(p+q, 'p_plus_q')
test_vector_point(p-q, 'p_minus_q')
test_vector_point(k*p, 'k_times_p')
test_vector_point(k1*p, 'k1_times_p')
test_vector_point(k2*p, 'k2_times_p')
test_vector_point(k1*p + k2*q, 'k1_times_p_plus_k2_times_q')
test_vector_scalar(k, 'k')
test_vector_scalar(k1, 'k1')
test_vector_scalar(k2, 'k2')

# small order point
test_vector_point((0, Fp(-1)), "p_order_2_1", False)
test_vector_point((0, 1, 0), "p_order_2_2", False, True)
test_vector_point((1, 0, 0), "p_order_2_3", False, True)

# # small x point
# xx = Fp(1)
# yy = sqrt((1-a_ed*xx**2)/(1-d_ed*xx**2))
# xx_ws = 29473314690983384253693538305232603171909485083805813020990978841153575434625
# yy_ws = 44104446522130670290417281820901020888519367242619929491530721381946199062276
# small_p_ws = E(xx_ws, yy_ws)
# assert small_p_ws.order() != r
# test_vector_point((xx, yy), "small_p", False)
# test_vector_point(2 * small_p_ws, "small_p_dbl")
