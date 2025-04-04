## Combining GLV and FakeGLV
We provide in `glv_fakeglv.py` a demonstration of a fast scalar multiplication computed using GLV and FakeGLV as presented in [this blogpost](https://ethresear.ch/t/fake-glv-you-dont-need-an-efficient-endomorphism-to-implement-glv-like-scalar-multiplication-in-snark-circuits/20394). 
While the computation is done in Python, it is aimed to be adapted for circuit integration. This example is conducted with Bandersnatch, but can be applied to any curve having a fast endomorphism (a-la-GLV), for example BLS12-381, BN254, or the Bitcoin's curve (depending on the context where an in-circuit verification is needed).

### How to use
```bash
cd ..
python -m example.glv_fakeglv --k 12345623456765432345676543234567876543456765434567865433456765433456765433456 --l 13456384938498474896734986734986734986739863749863746020720946720694720496702
```
### Scalar decomposition
This computes an example of scalar decomposition in dimension 4 for a random scalar $k$:
```python
Decomposition of k:
k = 12345623456765432345676543234567876543456765434567865433456765433456765433456 (253 bits)
u1 = 2232876471539877833 (61 bits)
u2 = 2399893344234450747 (62 bits)
v1 = -790970025150257770 (60 bits)
v2 = -9814722368297079247 (64 bits)
```
We provide an implementation that is not SageMath-dependent. Our implementation of LLL is very slow compared with SageMath (that wraps a C implementation), but in a ZK context, this computation is part of the witness and not computed in-circuit.

### Test of 4MSM
We verify that $[u_1]P + [u_2]\phi(P) - [v_1]Q - [v_2]\phi(Q) = 0$:
```
TEST of [u1]P + [u2]φ(P) - [v1]Q - [v2]φ(Q) == 0?		OK
```
We provide a quick benchmark (although this implementation is not optimized):
```
BENCHMARK of the methods:
GLV:       	    4.31 ms per iteration.
GLV+FakeGLV:	2.40 ms per iteration (44% faster).

```
We obtain a 44% improvement as expected from the theory: MSM(4,64) vs MSM(2,128).

### Test of 2MSM using GLV and FakeGLV
The code also provides a decomposition of

$$[k]P+[l]Q = R \iff [u_1]P + [u_2] \phi(P) + [v_1]Q + [v_2]\phi(Q) - [w_1]R -[w_2]\phi(R) = 0$$
```
Simultaneous decomposition of k and l:
k = 12345623456765432345676543234567876543456765434567865433456765433456765433456 (253 bits)
l = 347416144716927276873051607940243527430707856423175768934543413965032823901 (253 bits)
u1 = -8727058181421248280755369 (83 bits)
u2 = 598448052775363456475082 (79 bits)
v1 = 3688021787012064781328580 (82 bits)
v2 = 12046754130873887373428165 (84 bits)
w1 = -5038805389123913550157142 (83 bits)
w2 = 1209496161631686713732247 (81 bits)
```
This simultaneous decomposition is done using a lattice of dimension 5. We did not investigate further the implementation using a 6MSM yet, but we expect having a slight improvement of around 9% (MSM(6, 86) vs MSM(4, 128)).
