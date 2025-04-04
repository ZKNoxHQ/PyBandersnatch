# Implementation of GLV + Fake GLV
# find small u1,u2,v1,v2 such that s = (u1+λu2)/(v1+λv2) mod r
from time import time
from src.curve.edwards import Edwards
from src.field import Field
from example.lll import slow_reduction
import argparse


def splitGLVFakeGLV(k, λ, r, δ):
    """
    Decompose k = (u_1 + λ u_2) / (v_1 + v_2 λ) mod r using LLL with parameter δ
    Expected size |u_i|, |v_i| ≃ (r/α³)^(1/4)
    """

    basis = [
        [r, 0, 0, 0],
        [-λ, 1, 0, 0],
        [k, 0, 1, 0],
        [0, 0, -λ, 1],
    ]
    M = slow_reduction(basis, δ)
    return M[0]


def split2MSMGLVFakeGLV(k, l, λ, r, δ):
    """
    Decompose k = (u_1 + λ u_2) / (w_1 + λ w_2)
              l = (v_1 + λ v_2) / (w_1 + λ w_2)
    simultaneously using LLL with parameter δ
    Expected size |u_i|, |v_i|, |w_i| ≃ (r²λ/α⁴)**(1/5) ? Can be reduced ?
    """

    basis = [
        [r, 0, 0, 0, 0, 0],
        [0, 0, r, 0, 0, 0],
        [k, 0, l, 0, 1, 0],
        [-λ, 1, 0, 0, 0, 0],
        [0, 0, -λ, 1, 0, 0],
        [0, 0, 0, 0, -λ, 1],
    ]
    M = slow_reduction(basis, δ)
    return M[0]


def cli():
    parser = argparse.ArgumentParser(
        description="Demonstration of GLV+FakeGLV")
    # parser.add_argument("action", choices=[
    #                     "glv_fakeglv", "2msm_glv_fakeglv", "test_glv_fakeglv",], help="Action to perform")
    parser.add_argument("--k", type=int,
                        help="Provide a scalar to perform the scalar multiplication")
    parser.add_argument("--l", type=int,
                        help="Provide a second scalar to perform the multi scalar multiplication")

    args = parser.parse_args()

    if not args.k:
        print("Error: Provide --k")
        return

    # BANDERSNATCH ELLLIPTIC CURVE IN TWISTED EDWARDS MODEL
    F = Field(0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001)
    a = F(52435875175126190479447740508185965837690552500527637822603658699938581184508)
    d = F(45022363124591815672509500913686876175488063829319466900776701791074614335719)
    r = 0x1cfb69d4ca675f520cce760202687600ff8f87007419047174fd06b52876e7e1
    h = 4
    E = Edwards(a, d, r, h)
    λ = -8913659658109529928382530854484400854125314752504019737736543920008458395397

    # Verification of [k] P = Q
    P = E(F(0x1cc6ee38139c1c110223537a8ce79d067e58cc1067c6fbb7d8b3a1b08dfc8f08), F(
        0x70a5894a64445438d015ac32ba360f092cde44bab11fc2b7d4b5c0d216228cce), F(1))
    k = (args.k) % r
    Q = k*P

    # LLL parameter
    δ = 0.999
    α = (δ-1/4)
    theoretic_bound = round((r/α**3)**(1/4))

    (u1, u2, v1, v2) = splitGLVFakeGLV(k, λ, r, δ)
    assert (u1 + λ * u2 - k * v1 - k * λ * v2) % r == 0
    # in-circuit:
    # [u1]P + [u2]φ(P) - [v1] Q - [v2]φ(Q) == 0
    # where Q = [k]P is a hint.

    # Bound size
    assert max([abs(u1), abs(u2), abs(v1), abs(v2)]) < theoretic_bound

    print("Decomposition of k:")
    print("k = {} ({} bits)".format(k, len(bin(k))-2))
    print("u1 = {} ({} bits)".format(u1, len(bin(abs(u1)))-2))
    print("u2 = {} ({} bits)".format(u2, len(bin(abs(u2)))-2))
    print("v1 = {} ({} bits)".format(v1, len(bin(abs(v1)))-2))
    print("v2 = {} ({} bits)".format(v2, len(bin(abs(v2)))-2))

    print("\n\nTEST of [u1]P + [u2]φ(P) - [v1]Q - [v2]φ(Q) == 0?", end='\t\t')
    assert P.multi_scalar_mul_4(u1, P.φ(), u2, Q.neg(),
                                v1, Q.φ().neg(), v2) == E(0, 1, 1)
    if (P.multi_scalar_mul_4(u1, P.φ(), u2, Q.neg(),
                             v1, Q.φ().neg(), v2) == E(0, 1, 1)):
        print("OK")

    print("\n\nBENCHMARK of the methods:")
    niter = 500

    print("GLV:       ", end='\t')
    t = time()
    for i in range(niter):
        x = P.glv(k) + Q.neg()
    time_glv = (time() - t)/niter
    print("{:.2f} ms per iteration.".format(time_glv*1000))

    print("GLV+FakeGLV:", end='\t')
    t = time()
    for i in range(niter):
        x = P.multi_scalar_mul_4(u1, P.φ(), u2, Q.neg(), v1, Q.φ().neg(), v2)
    time_glv_fakeglv = (time()-t)/niter
    print("{:.2f} ms per iteration ({:.0f}% faster).".format(
        time_glv_fakeglv*1000, (time_glv - time_glv_fakeglv) / time_glv * 100))

    if args.l:

        # compute a 2MSM decomposition
        l = (args.l) % r
        (u1, u2, v1, v2, w1, w2) = split2MSMGLVFakeGLV(k, l, λ, r, δ)

        print("Simultaneous decomposition of k and l:")
        print("k = {} ({} bits)".format(k, len(bin(k))-2))
        print("l = {} ({} bits)".format(l, len(bin(k))-2))
        print("u1 = {} ({} bits)".format(u1, len(bin(abs(u1)))-2))
        print("u2 = {} ({} bits)".format(u2, len(bin(abs(u2)))-2))
        print("v1 = {} ({} bits)".format(v1, len(bin(abs(v1)))-2))
        print("v2 = {} ({} bits)".format(v2, len(bin(abs(v2)))-2))
        print("w1 = {} ({} bits)".format(w1, len(bin(abs(w1)))-2))
        print("w2 = {} ({} bits)".format(w2, len(bin(abs(w2)))-2))


cli()
