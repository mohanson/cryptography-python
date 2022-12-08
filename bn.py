P = 0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47


class Fp:
    def __init__(self, x):
        self.x = x % P

    def __repr__(self):
        return f'Fp(0x{self.x:064x})'

    def __eq__(self, other):
        return self.x == other.x

    def __add__(self, other):
        return Fp((self.x + other.x) % P)

    def __sub__(self, other):
        return Fp((self.x - other.x) % P)

    def __mul__(self, other):
        return Fp((self.x * other.x) % P)

    def __truediv__(self, other):
        return self * other ** -1

    def __pow__(self, other):
        return Fp(pow(self.x, other, P))

    def __neg__(self):
        return Fp(P - self.x)


def polyadd(c1, c2):
    return [c1[i] + c2[i] for i in range(len(c1))]


def polysub(c1, c2):
    return [c1[i] - c2[i] for i in range(len(c1))]


def polymul(c1, c2):
    p = [Fp(0) for _ in range(len(c1) + len(c2) - 1)]
    for i in range(len(c1)):
        for j in range(len(c2)):
            p[i+j] = p[i+j] + c1[i] * c2[j]
    return p


def polydivmod(c1, c2):
    # Algorithm: https://en.wikipedia.org/wiki/Polynomial_long_division
    # The code implementation is inspired by numpy.polynomial.polynomial.polydiv
    lc1 = len(c1)
    lc2 = len(c2)
    if lc1 < lc2:
        return [], c1
    if lc2 == 1:
        return [e / c2[0] for e in c1], []
    dif = lc1 - lc2
    scl = c2[-1]
    nc1 = c1.copy()
    nc2 = [e/scl for e in c2[:-1]]
    i = dif
    j = lc1 - 1
    while i >= 0:
        for k in range(lc2 - 1):
            nc1[i+k] -= nc2[k]*nc1[j]
        i -= 1
        j -= 1
    return [e/scl for e in nc1[j+1:]], nc1[:j+1]


def polydiv(c1, c2):
    return polydivmod(c1, c2)[0]


def polymod(c1, c2):
    return polydivmod(c1, c2)[1]


class Fp2:
    def __init__(self, coeffs):
        self.coeffs = coeffs
        self.mod = [Fp(1), Fp(0), Fp(1)]

    def __repr__(self):
        return f'Fp2({self.coeffs})'

    def __eq__(self, other):
        return self.coeffs == other.coeffs

    def __add__(self, other):
        return Fp2(polyadd(self.coeffs, other.coeffs))

    def __sub__(self, other):
        return Fp2(polysub(self.coeffs, other.coeffs))

    def __mul__(self, other):
        return Fp2(polymod(polymul(self.coeffs, other.coeffs), self.mod))

    def __truediv__(self, other):
        return Fp2(polydiv(self.coeffs, other.coeffs))


a = Fp2([Fp(3), Fp(0)])
b = Fp2([Fp(9), Fp(1)])
assert a + b == Fp2([Fp(12), Fp(1)])
assert a - b == Fp2([Fp(0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd41),
                     Fp(0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd46)])
assert a * b == Fp2([Fp(27), Fp(3)])
assert a / b == Fp2([Fp(0x2b149d40ceb8aaae81be18991be06ac3b5b4c5e559dbefa33267e6dc24a138e5),
                     Fp(0x009713b03af0fed4cd2cafadeed8fdf4a74fa084e52d1852e4a2bd0685c315d2)])

# def polyinv(a, p):
#     # https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Simple algebraic field extensions
#     t = [Fp(0) for _ in a];     newt = [Fp(1)] + [Fp(0) for _ in a[:-1]]
#     r = p.copy();  newr = a

#     while newr != [Fp(0) for _ in a]:
#         quotient = polydiv(r, newr)
#         r, newr = newr, polysub(r, polymul(quotient, newr))
#         t, newt = newt, polysub(t, polymul(quotient, newt))
#     # if degree(r) > 0:
#     #     return "Either p is not irreducible or a is a multiple of p"
#     return polymul(polydiv([Fp(1)] + [Fp(0) for _ in a[:-1]], r), t)

# print(polyinv(b.coeffs, [Fp(1), Fp(0), Fp(1)]))
