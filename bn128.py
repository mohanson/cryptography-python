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


class Fp12:
    def __init__(self, coeffs):
        self.coeffs = coeffs
        self.mod = [Fp(e) for e in [82, 0, 0, 0, 0, 0, -18, 0, 0, 0, 0, 0, 1]]  # w¹² - 18w⁶ + 82 = 0

    def __repr__(self):
        return f'Fp12({self.coeffs})'

    def __eq__(self, other):
        return self.coeffs == other.coeffs

    def __add__(self, other):
        return Fp12(polyadd(self.coeffs, other.coeffs))

    def __sub__(self, other):
        return Fp12(polysub(self.coeffs, other.coeffs))

    def __mul__(self, other):
        return Fp12(polymod(polymul(self.coeffs, other.coeffs), self.mod))

    def __truediv__(self, other):
        return self * Fp12(polyinv(other.coeffs, self.mod))

    def __pow__(self, other):
        if other == 0:
            return Fp12([Fp(e) for e in [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
        elif other == 1:
            return Fp12(self.coeffs)
        elif other % 2 == 0:
            return (self * self) ** (other // 2)
        else:
            return ((self * self) ** int(other // 2)) * self

    def __neg__(self):
        return Fp12([-c for c in self.coeffs])


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

a = Fp12([Fp(e) for e in [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4]])
b = Fp12([Fp(e) for e in [1, 2, 3, 4, 2, 2, 3, 4, 3, 2, 3, 4]])
assert a * b == Fp12([Fp(0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87bcee8),
                      Fp(0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87c404f),
                      Fp(0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87c5483),
                      Fp(0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87c6253),
                      Fp(0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87c9029),
                      Fp(0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87ceb79),
                      Fp(0x0000000000000000000000000000000000000000000000000000000000003327),
                      Fp(0x00000000000000000000000000000000000000000000000000000000000020b0),
                      Fp(0x0000000000000000000000000000000000000000000000000000000000001d44),
                      Fp(0x0000000000000000000000000000000000000000000000000000000000001a9a),
                      Fp(0x0000000000000000000000000000000000000000000000000000000000001321),
                      Fp(0x0000000000000000000000000000000000000000000000000000000000000438)])
assert a / b == Fp12([Fp(0x2cd098d90a44e85af5a59c764154b991c4afaec7bcd14a447961b0bd1d7d1d81),
                      Fp(0x2c9165696cd217b8d47c0e40cfe8f28a4c4cadccd8be630b036962ac44de2835),
                      Fp(0x000818c49d384098612bfbea1300d2fd2bffb545c7801dfd83422cd0baff6b6b),
                      Fp(0x16d6c80fddaa5d32be3c37ecb46ac06491cc2b10922e5b6bd929de6d21a15454),
                      Fp(0x2a1ec289458c31cdd73c149d0ba84dc40fe7a19e56c4d8216dad2cb56e1b616b),
                      Fp(0x0a1f284e8d2371a623abe8d6ac9a726239323f4fa9fce678670391061192afa7),
                      Fp(0x03c6a4a7b63aaa1e2d796e72b84ae8ffac2f00f53bf70a281703f900291604d8),
                      Fp(0x23d2fd15151a1b2d00962d62ef3ef7a6d2075a3cba98da69f808f90ebce3069d),
                      Fp(0x05012062c1750dfea579794c4a28d4dca8234045fe11e51363be1fa69fe01372),
                      Fp(0x0b3cab55aa0caf777904b997278f2034422b353d1be0795933a1d010c6b7adde),
                      Fp(0x05f55060cf8c3b66984c4212c3c5b693d61e6dbf46c46d92cf428317d91c007b),
                      Fp(0x2e0262f661d0c450b646bf2271134dbe8c306a6064e80aa15429af1a06d9a060)])
