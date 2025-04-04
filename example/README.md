## Combining GLV and FakeGLV
We provide in `glv_fakeglv.py` a demonstration of a fast scalar multiplication computed using GLV and FakeGLV. 
While the computations are done in Python, it is aimed to be adapted for circuit integration.

### How to use
```bash
cd ..
python -m example.glv_fakeglv
```
### Scalar decomposition
This computes an example of scalar decomposition in dimension 4 for a random scalar $k$:
```python
k = 8809196524735054409598625807987834789941239467291111440141961710399690321154 (253 bits)
u1 = -4721629758273561887 (64 bits)
u2 = 4445070398100683295 (62 bits)
v1 = -968749169646434063 (61 bits)
v2 = 2866665739561707568 (62 bits)
```
We provide an implementation that is not SageMath-dependent. Our implementation of LLL is very slow compared with SageMath (that wraps a C implementation), but in a ZK context, this computation is part of the witness and not computed in-circuit.

### Test of 4MSM
We verify that $[u_1]P + [u_2]φ(P) - [v_1]Q - [v_2]φ(Q) = 0$:
```
TEST of [u1]P + [u2]φ(P) - [v1]Q - [v2]φ(Q) == 0?		OK
```
We provide a quick benchmark (although this implementation is not optimized):
```
BENCHMARK of GLV+FakeGLV vs GLV:
GLV:       	    4.42 ms per iteration.
GLV+FakeGLV:	2.47 ms per iteration.
```
We obtain a 44% improvement as expected from the theory: MSM(4,64) vs MSM(2,128).