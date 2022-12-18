import bn128
import random
import polynomial

# [1] Polynomial commits: https://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf
# [2] https://www.youtube.com/watch?v=n4eiiCDhTes
# [3] https://www.youtube.com/watch?v=NVvNHe_RGZ8
# [4] https://copper-witch-857.notion.site/Polynomial-KZG-or-Kate-Commitment-DappLearning-Notes-fc426c8cb9a14878840852506865f13b


Fq = bn128.Fq
Fr = bn128.Fr
P1 = bn128.P1
G1 = bn128.G1
I1 = bn128.I1
G2 = bn128.G2
pairing = bn128.pairing

sk = Fr(random.randint(0, Fr.p - 1))
pk_g1 = [G1 * (sk**i) for i in range(10)]
pk_g2 = [G2 * (sk**i) for i in range(10)]

ax = [Fr(e) for e in [1,  2,  3,  4]]
ay = [Fr(e) for e in [4, 15, 40, 85]]
coeffs = polynomial.interp(ax, ay)
commit = I1
for i in range(len(coeffs)):
    commit = commit + pk_g1[i] * coeffs[i]

px = Fr(1)
py = Fr(4)
coeffs = polynomial.div(polynomial.sub(coeffs, [py]), [-px, Fr(1)])
proofs = I1
for i in range(len(coeffs)):
    proofs = proofs + pk_g1[i] * coeffs[i]

lhs = pairing(G2, commit - G1 * py)
rhs = pairing(pk_g2[1] - G2 * px, proofs)
assert lhs == rhs
