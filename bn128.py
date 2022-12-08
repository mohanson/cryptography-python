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


def polyclr(c1):
    for i in range(len(c1) - 1, -1, -1):
        if c1[i] != Fp(0):
            break
    return c1[:i+1]


def polydeg(c1):
    d = len(c1) - 1
    while c1[d] == Fp(0) and d:
        d -= 1
    return d


def polyadd(c1, c2):
    p = [Fp(0) for _ in range(max(len(c1), len(c2)))]
    for i, e in enumerate(c1):
        p[i] += e
    for i, e in enumerate(c2):
        p[i] += e
    return polyclr(p)


def polysub(c1, c2):
    p = [Fp(0) for _ in range(max(len(c1), len(c2)))]
    for i, e in enumerate(c1):
        p[i] += e
    for i, e in enumerate(c2):
        p[i] -= e
    return polyclr(p)


def polymul(c1, c2):
    p = [Fp(0) for _ in range(len(c1) + len(c2) - 1)]
    for i in range(len(c1)):
        for j in range(len(c2)):
            p[i+j] += c1[i] * c2[j]
    return polyclr(p)


def polydivmod(c1, c2):
    # Algorithm: https://en.wikipedia.org/wiki/Polynomial_long_division
    # The code implementation is inspired by numpy.polynomial.polynomial.polydiv
    lc1 = len(c1)
    lc2 = len(c2)
    if c2[-1] == Fp(0):
        raise ZeroDivisionError()
    if lc1 < lc2:
        return [Fp(0)], c1
    if lc2 == 1:
        return [e / c2[0] for e in c1], [Fp(0)]
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
    return [e/scl for e in nc1[j+1:]], polyclr(nc1[:j+1])


def polydiv(c1, c2):
    return polydivmod(c1, c2)[0]


def polymod(c1, c2):
    return polydivmod(c1, c2)[1]


def polyinv(c1, c2):
    newt, t = [Fp(1)], [Fp(0)]
    newr, r = c1, c2
    while polydeg(newr):
        quotient = polydiv(r, newr)
        r, newr = newr, polysub(r, polymul(newr, quotient))
        t, newt = newt, polysub(t, polymul(newt, quotient))
    return [e/newr[0] for e in newt[:polydeg(c2)]]


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
        return self * Fp2(polyinv(other.coeffs, self.mod))

    def __pow__(self, other):
        if other == 0:
            return Fp2([Fp(1), Fp(0)])
        elif other == 1:
            return Fp2(self.coeffs)
        elif other % 2 == 0:
            return (self * self) ** (other // 2)
        else:
            return ((self * self) ** int(other // 2)) * self

    def __neg__(self):
        return Fp2([-c for c in self.coeffs])


a = Fp2([Fp(3), Fp(0)])
b = Fp2([Fp(9), Fp(1)])
assert a + b == Fp2([Fp(12), Fp(1)])
assert a - b == Fp2([Fp(0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd41),
                     Fp(0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd46)])
assert a * b == Fp2([Fp(27), Fp(3)])
assert a / b == Fp2([Fp(0x2b149d40ceb8aaae81be18991be06ac3b5b4c5e559dbefa33267e6dc24a138e5),
                     Fp(0x009713b03af0fed4cd2cafadeed8fdf4a74fa084e52d1852e4a2bd0685c315d2)])
assert b ** 42 == Fp2([Fp(0x30644e72e131a029b85045b68181585aa4fe1cf011e7db175f2a2d79e17cfd47),
                       Fp(0x30644e72e131a029b85045b6818158302983bf245a26e5c61975a0d46d9cfd47)])
