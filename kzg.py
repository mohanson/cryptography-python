import bn128
import lagrange

lagrange.Fq = bn128.Fq
Fq = lagrange.Fq
Fr = bn128.Fr
lagrange.Foly.nil = Fr(0)
lagrange.Foly.one = Fr(1)

# God knowns secret
secret = Fr(123456)

# Create commitment
x = [Fr(e) for e in [1,  2,  3,  4]]
y = [Fr(e) for e in [4, 15, 40, 85]]
coeffs = lagrange.lagrange_foly(x, y)


def f(x):
    s = Fr(0)
    for i in range(len(coeffs)):
        s += (x ** i) * coeffs[i]
    return s


secvec = [bn128.G1 * (secret**i) for i in range(len(x))]
commit = bn128.Point(Fq(0), Fq(0))
for i in range(len(secvec)):
    commit += (secvec[i] * coeffs[i])
print('commit:', commit)

# Prover x=1 y=4
# f(x) - yᵢ = (x - xᵢ)q(x)
# q(x) = (f(x) - yᵢ)) / (x - xᵢ)


class FrPc(bn128.Pc):
    nil = Fr(0)
    one = Fr(1)


index = 0
coeffs_1 = [e for e in coeffs]
coeffs_1[0] = coeffs_1[0] - y[index]
coeffs_2 = [-x[index], Fr(1)]
qx_coeffs = FrPc.div(coeffs_1, coeffs_2)
print('qx_coeffs', qx_coeffs)
qs = Fr(0)
for i in range(len(qx_coeffs)):
    qs += qx_coeffs[i] * (secret ** i)
print('qs', qs)
prove = bn128.G1 * qs
print('prove', prove)

# Verifier e(C - yᵢ⋅G1, G2) = e(q(s)⋅G1, (s - xᵢ)⋅G2)
lhs = bn128.pairing(bn128.G2, commit - bn128.G1 * y[index])
rhs = bn128.pairing(bn128.G2 * secret - bn128.G2 * x[index], prove)
print(lhs)
print(rhs)
print(lhs == rhs)
