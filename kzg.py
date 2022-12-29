import bn128
import random
import polynomial

# [1] Polynomial commits: https://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf
# [2] https://www.youtube.com/watch?v=n4eiiCDhTes
# [3] https://www.youtube.com/watch?v=NVvNHe_RGZ8
# [4] https://copper-witch-857.notion.site/Polynomial-KZG-or-Kate-Commitment-DappLearning-Notes-fc426c8cb9a14878840852506865f13b
# [5] https://foresightnews.pro/article/detail/17988
# [6] https://dankradfeist.de/ethereum/2021/10/13/kate-polynomial-commitments-mandarin.html


Fr = bn128.Fr
G1 = bn128.G1
I1 = bn128.I1
I2 = bn128.I2
G2 = bn128.G2
pairing = bn128.pairing

sk = Fr(random.randint(0, Fr.p - 1))
pk_g1 = [G1 * (sk**i) for i in range(10)]
pk_g2 = [G2 * (sk**i) for i in range(10)]

ax = [Fr(e) for e in [0,  1,  2,  3]]
ay = [Fr(e) for e in [4, 15, 40, 85]]
fx = polynomial.lagrange(ax, ay)
commit = sum([pk_g1[i] * fx[i] for i in range(len(fx))], start=I1)

# 单个证明
px = Fr(1)
py = Fr(15)
qx = polynomial.div(polynomial.sub(fx, [py]), [-px, Fr(1)])
proofs = sum([pk_g1[i] * qx[i] for i in range(len(qx))], start=I1)
lhs = pairing(G2, commit - G1 * py)
rhs = pairing(pk_g2[1] - G2 * px, proofs)
assert lhs == rhs

# 批量证明
px = [Fr(e) for e in [0,  1]]
py = [Fr(e) for e in [4, 15]]
ix = polynomial.lagrange(px, py)
zx = polynomial.zerofier(px)
qx = polynomial.div(polynomial.sub(fx, ix), zx)
proofs = sum([pk_g1[i] * qx[i] for i in range(len(qx))], start=I1)
ix_sg1 = sum([pk_g1[i] * ix[i] for i in range(len(ix))], start=I1)
zx_sg2 = sum([pk_g2[i] * zx[i] for i in range(len(zx))], start=I2)
lhs = pairing(G2, commit - ix_sg1)
rhs = pairing(zx_sg2, proofs)
assert lhs == rhs
