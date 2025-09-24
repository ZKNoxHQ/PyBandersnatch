from hashlib import sha512


Fp = GF((1 << 255)-19)
Fpxyz = Fp['x', 'y', 'z']
Fpxyz = FractionField(Fpxyz)
x, y, z = Fpxyz.gens()
r = 7237005577332262213973186563042994240857116359379907606001950938285454250989
h = 4

# Ed25519, Twisted-Edwards model
a_ed, d_ed = Fp(-1), Fp(-121665/121666)


# Conversion Montgomery <-> (twisted) Edwards
sqrt_m486664 = Fp(
    51042569399160536130206135233146329284152202253034631822681833788666877215207)
assert sqrt_m486664**2 + 486664 == 0
ed2mg_x = (1+y) / (1-y)
ed2mg_y = sqrt_m486664 * (1+y)/(1-y) / x
mg2ed_x = sqrt_m486664 * x/y
mg2ed_y = (x-1)/(x+1)

a_mg, b_mg = Fp(486662), Fp(1)  # 2*(a_ed+d_ed)/(a_ed-d_ed), 4/(a_ed-d_ed)

# Conversion Montgomery <-> Weierstrass
mg2ws_x = x/b_mg + a_mg / (3*b_mg)
mg2ws_y = y/b_mg
ws2mg_x = b_mg * (x-(a_mg/(3*b_mg)))
ws2mg_y = b_mg * y
a_ws, b_ws = (3-a_mg**2)/(3*b_mg**2), (2*a_mg**3 - 9 * a_mg) / (27*b_mg**3)


# Taken from RFC-8032
x_ed, y_ed = 15112221349535400772501151409588531511454012693041857206046113283949847762202, 46316835694926478169428394003475163141307993866256225615783033603165251855960
assert a_ed * x_ed**2 + y_ed**2 == 1 + d_ed * x_ed**2 * y_ed**2
x_mg, y_mg = ed2mg_x(y=y_ed), ed2mg_y(x=x_ed, y=y_ed)
assert b_mg * y_mg**2 == x_mg**3 + a_mg * x_mg**2 + x_mg
assert mg2ed_x(x=x_mg, y=y_mg) == x_ed
assert mg2ed_y(x=x_mg, y=y_mg) == y_ed
x_ws, y_ws = mg2ws_x(x=x_mg, y=y_mg), mg2ws_y(x=x_mg, y=y_mg)
assert y_ws**2 == x_ws**3 + a_ws*x_ws + b_ws
assert ws2mg_x(x=x_ws, y=y_ws) == x_mg
assert ws2mg_y(x=x_ws, y=y_ws) == y_mg

E = EllipticCurve([a_ws, b_ws])
p = E(x_ws, y_ws)
assert p.order() == r

# Scalar from RFC 8032
sca = bytes.fromhex(
    "9d61b19deffd5a60ba844af492ec2cc44449c5697b326919703bac031cae7f60")
pk = bytes.fromhex(
    "d75a980182b10ab7d54bfed3c964073a0ee172f3daa62325af021a68f707511a")


def secret_expand(secret):
    if len(secret) != 32:
        raise Exception("Bad size of private key")
    h = sha512(secret).digest()
    a = int.from_bytes(h[:32], "little")
    a &= (1 << 254) - 8
    a |= (1 << 254)
    return (a, h[32:])


(sk, _) = secret_expand(sca)

q = sk * p

q_x, q_y = q[0]/q[2], q[1]/q[2]
q_x_mg = ws2mg_x(x=q_x, y=q_y)
q_y_mg = ws2mg_y(x=q_x, y=q_y)
assert b_mg * q_y_mg**2 == q_x_mg**3 + a_mg * q_x_mg**2 + q_x_mg

q_x_ed = mg2ed_x(x=q_x_mg, y=q_y_mg)
q_y_ed = mg2ed_y(x=q_x_mg, y=q_y_mg)
assert a_ed * q_x_ed**2 + q_y_ed**2 == 1 + d_ed * q_x_ed**2 * q_y_ed**2

recover_pk = (ZZ(q_y_ed) | ((ZZ(q_x_ed) & 1) << 255)).to_bytes(32, 'little')
assert recover_pk == pk
