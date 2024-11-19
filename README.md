# Bxndxrsnxtch
A constant time x-only implementation of Bandersnatch elliptic curve.

## Context
[Bandersnatch](https://eprint.iacr.org/2021/1152.pdf) is an elliptic curve designed for zero-knowledge proof computations.
It allows efficient scalar multiplications using the [GLV method](https://www.iacr.org/archive/crypto2001/21390189.pdf).
We implement the XZ-only Montgomery representation of the curve in order to improve the scalar multiplication efficiency.
The project is written in Python.
The integer arithmetic is computed using `gmpy2`, a wrapper to the efficient `gmp` library written in `C`.

## Test
To generate the test vectors:
```
sage field_test_vectors.sage > field_test_vectors.py
sage curve_test_vectors.sage > curve_test_vectors.py
```

Tests for `Field` and `Curve`can be computed using:
```
python main.py --test
```
Benchmarks for GLV vs scalar multiplication can be computed using:
```
python main.py --bench
```