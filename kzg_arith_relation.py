import bn128
import hashlib
import random
import polynomial

Fr = bn128.Fr
G1 = bn128.G1
I1 = bn128.I1
I2 = bn128.I2
G2 = bn128.G2
pairing = bn128.pairing


def create_witness(pk, fx, x):
    y = polynomial.evaluate(fx, x)
    qx = polynomial.div(polynomial.sub(fx, [y]), [-x, Fr(1)])
    proofs = sum([pk[i] * qx[i] for i in range(len(qx))], start=I1)
    return x, y, proofs


sk = Fr(random.randint(0, Fr.p - 1))
pk_g1 = [G1 * (sk**i) for i in range(10)]
pk_g2 = [G2 * (sk**i) for i in range(10)]

# Prove that f(x) + g(x) = h(x)

fx = [Fr(e) for e in [1, 2, 3]]
gx = [Fr(e) for e in [4, 5, 6]]
hx = [Fr(e) for e in [5, 7, 9]]
fx_commit = sum([pk_g1[i] * fx[i] for i in range(len(fx))], start=I1)
gx_commit = sum([pk_g1[i] * gx[i] for i in range(len(gx))], start=I1)
hx_commit = sum([pk_g1[i] * hx[i] for i in range(len(hx))], start=I1)
m = hashlib.sha256()
m.update(fx_commit.x.x.to_bytes(32, 'little'))
m.update(fx_commit.y.x.to_bytes(32, 'little'))
m.update(gx_commit.x.x.to_bytes(32, 'little'))
m.update(gx_commit.y.x.to_bytes(32, 'little'))
m.update(hx_commit.x.x.to_bytes(32, 'little'))
m.update(hx_commit.y.x.to_bytes(32, 'little'))
z = Fr(int.from_bytes(m.digest(), 'little'))
_, fz, fx_proofs = create_witness(pk_g1, fx, z)
_, gz, gx_proofs = create_witness(pk_g1, gx, z)
_, hz, hx_proofs = create_witness(pk_g1, hx, z)
assert fz + gz == hz
assert pairing(G2, fx_commit - G1 * fz) == pairing(pk_g2[1] - G2 * z, fx_proofs)
assert pairing(G2, gx_commit - G1 * gz) == pairing(pk_g2[1] - G2 * z, gx_proofs)
assert pairing(G2, hx_commit - G1 * hz) == pairing(pk_g2[1] - G2 * z, hx_proofs)
