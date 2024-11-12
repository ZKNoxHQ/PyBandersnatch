# FaultyECDSA
An implementation of ECDSA allowing fault injection.

## Context
[Bandersnatch](https://eprint.iacr.org/2021/1152.pdf) is an elliptic curve designed for zero-knowledge proof computations.
It allows efficient scalar multiplications using the [GLV method](https://www.iacr.org/archive/crypto2001/21390189.pdf).
We implement the XZ-only Montgomery representation of the curve in order to improve the scalar multiplication efficiency.
The project is written in Python.
The integer arithmetic is computed using `gmpy2`, a wrapper to the efficient `gmp` library written in `C`.

## ECDSA efficiency
We implement the elliptic curve digital signal algorithm (ref?). We provide some benchmarks TODO.

## Fault injection
We provide a demonstration of the fault injection attack in this context. TODO.
This can be done using the command
```
python SOMETHING
```

## Test
Many test functions are written and can be computed using:
```
python -m unittest field_test curve_test -v
```