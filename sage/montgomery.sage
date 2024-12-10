ENDOMORPHISM = 1
if ENDOMORPHISM:
    # Bandersnatch, Montgomery model
    load("sage/φ.sage")
    φ_coefficients = montgomery_φ_coefficients
    Fp = GF(
        52435875175126190479447740508185965837690552500527637822603658699938581184513)
    A = Fp(0x4247698f4e32ad45a293959b4ca17afa4a2d2317e4c6ce5023e1fd63d1b5de98)
    B = Fp(5)
    r = 13108968793781547619861935127046491459309155893440570251786403306729687672801
    h = 4
else:
    # Curve25519, Montgomery model
    # [α,β,γ, ????????] defining φ = Id
    φ_coefficients = [1, 0, 0, 0, 0, 0, 0, 1]
    Fp = GF((1 << 255)-19)
    A = Fp(486662)
    B = Fp(1)
    r = (1 << 252)+27742317777372353535851937790883648493
    h = 8
    glv_params = [-1, 1, (r-1)//2]

# Weierstrass model
a = B**2 * (1-A**2/3)
b = B**3 * A * (2*A**2/9 - 1) / 3
E = EllipticCurve(Fp, [a, b])
# assert (E.j_invariant() == 8000)
assert (E.order() % r == 0)


def to_mg(P):
    if P[2] == 0:
        return 1, 0
    u = P[0]/P[2]
    v = P[1]/P[2]
    x = u/B - A/3
    y = v/B**2
    return x, y


if ENDOMORPHISM:
    # Bandersnatch
    p = E(31511963179209183026886029814959507395230513391536014203721350106469568871776,
          45347120062487836513813256222005391829075297965413648488198604153937949600247)
    q = E(20252373884274151187306374916054971403178400027240398097441693977206289492028,
          51776867226593987565156122032412653740536539297269422625045303422905294395870)
    λ = -8913659658109529928382530854484400854125314752504019737736543920008458395397
else:
    # Curve25519
    p = E(50985278150407432344796787702948612699476491671352921418144363141586907039370,
          15253002477052688986620419946290387236033221949375887382230120676776923013400)
    q = E(15681357382615180190462123999934171085226024319595864033686213637078275740748,
          47407999444187020738561942407576495663998743082126352560608723174318508905182)
    λ = 1

# p = h*E.random_point()
# print("p = E({}, {})".format(p[0], p[1]))
# q = h*E.random_point()
# print("q = E({}, {})".format(q[0], q[1]))
# # λ is the eigenvalue of φ

k = 11997154529596648729624281997554038960651754640906483911385998427296165917073
k1 = -45894336995428141233269187797940484884
k2 = 8683555061824981937504960049179714114


def φ(p):
    if ENDOMORPHISM:
        # endomorphism in affine coordinates as in 2021/1152.pdf
        x_p, y_p = p[0]/p[2], p[1]/p[2]
        E = p.curve()
        alpha = E.division_polynomial(2).roots()[0][0]
        P = E.lift_x(alpha)
        phi0, phi1 = E.isogeny(P)
        E2 = E.isogeny_codomain(P)
        # Isomorphism
        u = (E.a4()/E2.a4()).sqrt().sqrt()
        assert u**4 == E.a4()/E2.a4() and u**6 == E.a6()/E2.a6()
        rX = phi0(y=1) * u**2
        sX = phi1(y=1) * u**3
        return E(rX(x=x_p, y=y_p), y_p * sX(x=x_p, y=y_p), 1)
    else:
        return p


def φ_minus_one(p):
    return φ(p)-p


# def φ_minus_one(p):
#     # see `φ.sage`
#     X, Y = to_mg(p)
#     α = Fp(13017314467421381532402061398313046228820690393386411611562176812113295071440)
#     β = Fp(14989411347484419666605643019079533103863186413725217032868654387860539633484)
#     γ = Fp(39953720565912266872856944794434720047230584117801669040511822283402326025498)
#     return α*X*(X+β)**2/(X+γ)**2


# def φ(p):
#     # see r_mg_XY and s_mg_XY in `φ.sage`
#     X, Y = to_mg(p)
#     a1 = 26217937587563095239723870254092982918845276250263818911301829349969290592256
#     a2 = 14989411347484419663140498193005880785086916883037474254598401919095177670475
#     a3 = 14989411347484419663140498193005880785086916883037474254598401919095177670477
#     return a1 * (X+a2) * (X+a3)/X

# GENERATION OF TEST VECTORS


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
print("from src.curve.montgomery import Montgomery")
print("F = Field({})".format(Fp.characteristic()))
print("a = F({})".format(A))
print("b = F({})".format(B))
print("r = {}".format(hex(r)))
print("h = {}".format(h))
print("φ_coefficients = {}".format(φ_coefficients))
print("glv_params = {}".format(glv_params))
print("E = Montgomery(a, b, r, h, φ_coefficients, glv_params)")
print("test_vectors = {}")
test_vector_point(p, 'p')
test_vector_point(q, 'q')
test_vector_point(2*p, 'p_double')
test_vector_point(p+q, 'p_plus_q')
test_vector_point(p-q, 'p_minus_q')
test_vector_point(φ(p), 'φ_p')
test_vector_point(φ_minus_one(p), 'φ_minus_one_p')
test_vector_point(k*p, 'k_times_p')
test_vector_point(k1*p, 'k1_times_p')
test_vector_point(k2*p, 'k2_times_p')
test_vector_point(k1*p + k2*q, 'k1_times_p_plus_k2_times_q')
test_vector_scalar(k, 'k')
test_vector_scalar(k1, 'k1')
test_vector_scalar(k2, 'k2')
test_vector_scalar(λ, 'λ')
