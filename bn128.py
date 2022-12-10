P = 0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47


class Fg:
    # Galois field. In mathematics, a finite field or Galois field is a field that contains a finite number of elements.
    # As with any field, a finite field is a set on which the operations of multiplication, addition, subtraction and
    # division are defined and satisfy certain basic rules.

    p = 0

    def __init__(self, x):
        self.x = x % self.p

    def __repr__(self):
        return f'Fp(0x{self.x:064x})'

    def __eq__(self, other):
        return self.x == other.x

    def __add__(self, other):
        return self.__class__((self.x + other.x) % self.p)

    def __sub__(self, other):
        return self.__class__((self.x - other.x) % self.p)

    def __mul__(self, other):
        return self.__class__((self.x * other.x) % self.p)

    def __truediv__(self, other):
        return self * other ** -1

    def __pow__(self, other):
        return self.__class__(pow(self.x, other, self.p))

    def __neg__(self):
        return self.__class__(self.p - self.x)


Fg.p = 23
assert Fg(12) + Fg(20) == Fg(9)
assert Fg(8) * Fg(9) == Fg(3)
assert Fg(8) ** -1 == Fg(3)


class Fp(Fg):
    p = P

    def __repr__(self):
        return f'Fp(0x{self.x:064x})'


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


def polyext(c1, sz):
    p = [Fp(0) for _ in range(sz)]
    for i, e in enumerate(c1):
        p[i] = e
    return p


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
    # Algorithm: https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
    newt, t = [Fp(1)], [Fp(0)]
    newr, r = c1, c2
    while polydeg(newr):
        quotient = polydiv(r, newr)
        r, newr = newr, polysub(r, polymul(newr, quotient))
        t, newt = newt, polysub(t, polymul(newt, quotient))
    return polyclr([e/newr[0] for e in newt[:polydeg(c2)]])


class Fpx:
    # A class for elements in polynomial extension fields

    degree = 0
    p = []

    def __init__(self, coeffs):
        assert len(coeffs) == self.degree
        self.coeffs = coeffs

    def __repr__(self):
        return f'Fpx({self.coeffs})'

    def __eq__(self, other):
        return self.coeffs == other.coeffs

    def __add__(self, other):
        return self.__class__(polyext(polyadd(self.coeffs, other.coeffs), self.degree))

    def __sub__(self, other):
        return self.__class__(polyext(polysub(self.coeffs, other.coeffs), self.degree))

    def __mul__(self, other):
        return self.__class__(polyext(polymod(polymul(self.coeffs, other.coeffs), self.p), self.degree))

    def __truediv__(self, other):
        return self * self.__class__(polyext(polyinv(other.coeffs, self.p), self.degree))

    def __pow__(self, other):
        if other == 0:
            return self.__class__([Fp(1)] + [Fp(0) for _ in range(self.degree - 1)])
        if other == 1:
            return self.__class__([e for e in self.coeffs])
        if other % 2 == 0:
            return (self * self) ** (other // 2)
        return (self * self) ** (other // 2) * self

    def __neg__(self):
        return self.__class__([-c for c in self.coeffs])


class Fp2(Fpx):
    degree = 2
    p = [Fp(e) for e in [1, 0, 1]]  # i² + 1 = 0


class Fp12(Fpx):
    degree = 12
    p = [Fp(e) for e in [82, 0, 0, 0, 0, 0, -18, 0, 0, 0, 0, 0, 1]]  # w¹² - 18w⁶ + 82 = 0


N = 21888242871839275222246405745257275088548364400416034343698204186575808495617

# Curve order should be prime
assert pow(2, N, N) == 2
# Curve order should be a factor of P**12 - 1
assert (P ** 12 - 1) % N == 0

# Curve is y**2 = x**3 + 3
B = Fp(3)
# Twisted curve over FQ**2
B2 = Fp2([Fp(3), Fp(0)]) / Fp2([Fp(9), Fp(1)])
# Extension curve over FQ**12; same b value as over FQ
B12 = Fp12([Fp(e) for e in [3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
