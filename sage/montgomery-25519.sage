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


def to_mg(P):
    u = P[0]/P[2]
    v = P[1]/P[2]
    x = u/B - A/3
    y = v/B**2
    return x, y


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


def test_vector_point(p, name, ws=True):
    if ws:
        [x, y] = to_mg(p)
    else:
        x = p
    print("test_vectors['{}'] = E(F({}), 1)".format(name, hex(x), 1))


def test_vector_scalar(k, name):
    print("test_vectors['{}'] = {}".format(name, hex(k)))


print("# File generated using `sage sage/montgomery.sage > tests/vectors/montgomery.py`.")
print("from src.field import Field")
print("from src.curve.montgomery25519 import Montgomery25519")
print("F = Field({})".format(Fp.characteristic()))
print("a = F({})".format(A))
print("b = F({})".format(B))
print("r = {}".format(hex(r)))
print("h = {}".format(h))
print("E = Montgomery25519(a, b, r, h)")
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
