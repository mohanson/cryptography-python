import bn128
import random
import polynomial

# Polynomial commits: https://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf

Fq = bn128.Fq
Fr = bn128.Fr
P1 = bn128.P1
G1 = bn128.G1
I1 = bn128.I1
G2 = bn128.G2
pairing = bn128.pairing

sk = Fr(random.randint(0, Fr.p - 1))
pk = [G1 * (sk**i) for i in range(10)]

ax = [Fr(e) for e in [1,  2,  3,  4]]
ay = [Fr(e) for e in [4, 15, 40, 85]]
coeffs = polynomial.interp(ax, ay)
commit = P1(Fq(0), Fq(0))
for i in range(len(coeffs)):
    commit = commit + pk[i] * coeffs[i]

px = Fr(1)
py = Fr(4)
coeffs_1 = [e for e in coeffs]
coeffs_1[0] = coeffs_1[0] - py
coeffs_2 = [-px, Fr(1)]
coeffs_q = polynomial.div(coeffs_1, coeffs_2)
proofs = I1
for i in range(len(coeffs_q)):
    proofs = proofs + pk[i] * coeffs_q[i]

lhs = pairing(G2, commit - G1 * py)
rhs = pairing(G2 * sk - G2 * px, proofs)
assert lhs == rhs
