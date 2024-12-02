# Endomorphism in Montgomery model

## Group law
We follow Sections 3.2, 3.3 of [this paper](https://eprint.iacr.org/2017/212.pdf) for the $XZ$-only arithmetic of Bandersnatch:

$$\begin{array}{rl}
X_{[2]P} &= (X_P+ Z_P)²(X_P-Z_P)²,\\
Z_{[2]P} &= (4X_PZ_P)((X_P-Z_P)² + ((A+2)/4 (4X_PZ_P))
\end{array}$$

$$\begin{array}{rl}
X_{P+Q} &= Z_{P-Q} [(X_P-Z_P)(X_Q+Z_Q) + (X_P+Z_P)(X_Q-Z_Q)]², \\
Z_{P+Q} &= X_{P-Q} [(X_P-Z_P)(X_Q+Z_Q) - (X_P+Z_P)(X_Q-Z_Q)]²
\end{array}$$

The choice of model influences the cost of additions and doubling of the curve. Several models are presented in [this website](https://hyperelliptic.org/EFD/g1p/index.html). In particular,
* The **twisted Edwards** model computes a doubling in **7** multiplications and an addition in **10** multiplications,
* The **short Weierstrass** model computes a doubling in **10** multiplications and an addition in **12** multiplications,
* The **Montgomery** model computes a doubling in **4** multiplications and an addition in **6** multiplications. This model is the most efficient, but the addition requires the knowledge of $x(P-Q)$ in order to compute $x(P+Q)$. Also, this model represents the points of the curve up to a sign.

## GLV in XZ-only coordinates
Bandersnatch is designed so that it has an efficient endomorphism $\psi$. In $XZ$-only coordinates, this endomorphism is given in [this paper](https://eprint.iacr.org/2021/1152.pdf) (page 6):

$$\begin{array}{rl}
X_{φ(P)} &= -(X_P-Z_P)² - (A+2)X_PZ_P,\\
Z_{φ(P)} &= 2X_PZ_P
\end{array}$$

In order to implement GLV, the scalar multiplication is computed by:
1. Decomposing $k = k_1+\sqrt{-2} k_2$ where $\sqrt{-2} \in\mathbb Z / r\mathbb Z$ where $r$ is the scalar field of Bandersnatch. The scalars $k_1$ and $k_2$ can be found half of the size of $k$.
2. Precompute $P$, $φ(P)$ and $P+φ(P)$,
3. Compute $[k]P = [k_1]P + [k_2] φ(P)$ using a double scalar multiplication.

The overall gain is that we compute $\log_2(k)/2$ doubling instead of $\log_2(k)$ in the naive multiplication case.

In the context of the Montgomery model, computing $P+φ(P)$ is not straightforward, as the addition requires the knowledge of $P-φ(P)$ in the addition formula above. The (non-differential) addition formula is provided in [this paper](https://www.iacr.org/archive/eurocrypt2014/84410275/84410275.pdf) (page 8):

$$X(P±Q) = \frac{B(X_PY_Q∓X_QY_P)²}{X_PX_Q(X_P-X_Q)²}.$$

Using $φ(P)$ instead of $Q$, we obtain:

$$X(P±φ(P)) = \frac{B(X_PY_{φ(P)}∓X_{φ(P)}Y_P)²}{X_PX_{φ(P)}(X_P-X_{φ(P)})²}.$$

As $φ$ is an endomorphism. the value of $X_{φ(P)}$ and $Y_{φ(P)}$ are given by

$$\begin{array}{rl}
X_{φ(P)} &= r(X_P),\\
Y_{φ(P)} &= Y_P\cdot r(X_P)
\end{array}$$

for a rational function $r(X_P) \in \mathbb F_p(X_P)$.
This function can be derived using SageMath (see [here](sage/φ.sage)).
From this, the addition formula simplifies when we add $P$ and $φ(P)$:

$$
X(P±φ(P)) = \frac{B(X_PY_P\cdot r(X_P)∓r(X_P)Y_P)²}{X_P\cdot r(X_P)(X_P-r(X_P))²}
= \frac{BY_P²\cdot r(X_P)(X_P∓1)²}{X_P(X_P-r(X_P))²}$$

Finally, $BY_P² = X_P³+AX_P²+X_P$ and we obtain a formula depending only on $X_P$:

$$X(P±φ(P)) = \frac{(X_P³+AX_P²+X_P)\cdot r(X_P)(X_P∓1)²}{X_P(X_P-r(X_P))²}$$

Replacing $r(X_P)$, we obtain a simple expression with three field elements $α,β,γ$:

$$X(P-φ(P)) = \frac{αX_P(X_P+Z_Pβ)²}{Z_P(X_P+γZ_P)²}$$

This expression can be precomputed with few multiplications in order to apply the GLV technique.
