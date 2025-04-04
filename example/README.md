## Combining GLV and FakeGLV
We provide in `glv_fakeglv.py` a demonstration of a fast scalar multiplication computed using GLV and FakeGLV. 
While the computations are done in Python, it is aimed to be adapted for circuit integration.

### How to use
```bash
cd ..
python -m example.glv_fakeglv --k 12345623456765432345676543234567876543456765434567865433456765433456765433456 
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
We verify that $[u_1]P + [u_2]φ(P) - [v_1]Q - [v_2]φ(Q) = 0$:
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