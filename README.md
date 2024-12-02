# PyBandersnatch
A Python implementation of xECDH and EdDSA with **[Bandersnatch](https://eprint.iacr.org/2021/1152.pdf)**:
* **xECDA**: implementation in (XZ-only) Montgomery model, following [RFC 7748](https://datatracker.ietf.org/doc/html/rfc7748) with optimized GLV scalar multiplication,
* **EdDSA**: implementation in (twisted) Edwards model, following [RFC 8032](https://datatracker.ietf.org/doc/html/rfc8032) with optimized GLV scalar multiplication.

This implementationd does not has `sage` dependencies. The integer arithmetic is computed using `gmpy2`, a wrapper to `gmp` written in `C`.

## Context
[Bandersnatch](https://eprint.iacr.org/2021/1152.pdf) is an elliptic curve designed for zero-knowledge proof computations.
It allows ZK proofs computed with BLS12-381, a pairing-friendly curve designed for ZK applications.

## How to use

### Install
```
make install
```
As of today, for mac OS, `gmpy2` is not supported by brew, so it is required to run a venv to run the library:

```
python3 -m venv path/to/venv 
source path/to/venv/bin/activate
python3 -m venv path/to/venv
pip3 install --global-option=build_ext --global-option="-I/opt/homebrew/include/" --global-option="-L/opt/homebrew/lib/" gmpy2
```

### Generate test vectors (optional, requires `sage`)
```
make gen_test_vec
```

### Tests
```
make test
```
For testing only one specific function (for example here `test_sign_verify`):
```
make test TEST=test_eddsa.TestEdDSA.test_sign_verify
```

### Benchmarks
```
make benchmark
```
