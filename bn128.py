import polynomial


class Fp:
    # Galois field. In mathematics, a finite field or Galois field is a field that contains a finite number of elements.
    # As with any field, a finite field is a set on which the operations of multiplication, addition, subtraction and
    # division are defined and satisfy certain basic rules.
    #
    # https://www.cs.miami.edu/home/burt/learning/Csc609.142/ecdsa-cert.pdf
    # Don Johnson, Alfred Menezes and Scott Vanstone, The Elliptic Curve Digital Signature Algorithm (ECDSA)
    # 3.1 The Finite Field Fp

    p = 0

    def __init__(self, x):
        self.x = x % self.p

    def __repr__(self):
        return f'Fp(0x{self.x:064x})'

    def __eq__(self, data):
        assert self.p == data.p
        return self.x == data.x

    def __add__(self, data):
        assert self.p == data.p
        return self.__class__((self.x + data.x) % self.p)

    def __sub__(self, data):
        assert self.p == data.p
        return self.__class__((self.x - data.x) % self.p)

    def __mul__(self, data):
        assert self.p == data.p
        return self.__class__((self.x * data.x) % self.p)

    def __truediv__(self, data):
        return self * data ** -1

    def __pow__(self, data):
        return self.__class__(pow(self.x, data, self.p))

    def __pos__(self):
        return self

    def __neg__(self):
        return self.__class__(self.p - self.x)

    @classmethod
    def nil(cls):
        return cls(0)

    @classmethod
    def one(cls):
        return cls(1)


if __name__ == '__main__':
    Fp.p = 23
    assert Fp(12) + Fp(20) == Fp(9)
    assert Fp(8) * Fp(9) == Fp(3)
    assert Fp(8) ** -1 == Fp(3)
    Fp.p = 0

# Prime of finite field.
P = 0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47
# The order n of G.
N = 0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001
assert pow(2, N, N) == 2
assert (P ** 12 - 1) % N == 0


class Fq(Fp):

    p = P

    def __repr__(self):
        return f'Fq(0x{self.x:064x})'


class Fr(Fp):

    p = N

    def __repr__(self):
        return f'Fr(0x{self.x:064x})'


if __name__ == '__main__':
    Fp.p = 13
    assert polynomial.interp([Fp(1), Fp(4)], [Fp(6), Fp(2)]) == [Fp(3), Fp(3)]
    Fp.p = 0


class Pa:
    a = None
    b = None
    i = None

    def __init__(self, x, y):
        if x != self.i[0] or y != self.i[1]:
            assert y ** 2 == x ** 3 + self.a * x + self.b
        self.x = x
        self.y = y

    def __repr__(self):
        return f'Pa({self.x}, {self.y})'

    def __eq__(self, data):
        return self.x == data.x and self.y == data.y

    def __add__(self, data):
        # https://www.cs.miami.edu/home/burt/learning/Csc609.142/ecdsa-cert.pdf
        # Don Johnson, Alfred Menezes and Scott Vanstone, The Elliptic Curve Digital Signature Algorithm (ECDSA)
        # 4.1 Elliptic Curves Over Fp
        if self.x == self.i[0] and self.y == self.i[1]:
            return data
        if data.x == self.i[0] and data.y == self.i[1]:
            return self
        if self.x == data.x and self.y == -data.y:
            return self.__class__(self.i[0], self.i[1])
        x1, x2 = self.x, data.x
        y1, y2 = self.y, data.y
        if self.y == data.y:
            s = (x1 * x1 + x1 * x1 + x1 * x1 + self.a) / (y1 + y1)
        else:
            s = (y2 - y1) / (x2 - x1)
        x3 = s * s - x1 - x2
        y3 = s * (x1 - x3) - y1
        return self.__class__(x3, y3)

    def __sub__(self, data):
        return self + data.__neg__()

    def __mul__(self, k):
        # Point multiplication: Double-and-add
        # https://en.wikipedia.org/wiki/Elliptic_curve_point_multiplication
        n = k.x
        result = self.__class__(self.i[0], self.i[1])
        addend = self
        while n:
            b = n & 1
            if b == 1:
                result += addend
            addend = addend + addend
            n = n >> 1
        return result

    def __pos__(self):
        return self

    def __neg__(self):
        return self.__class__(self.x, -self.y)


class P1(Pa):
    a = Fq(0)
    b = Fq(3)
    i = [
        Fq(0),
        Fq(0)
    ]


G1 = P1(Fq(1), Fq(2))
I1 = P1(Fq(0), Fq(0))

if __name__ == '__main__':
    assert G1 * Fr(2) + G1 + G1 == G1 * Fr(4)
    assert G1 + G1 != G1
    assert G1 * Fr(9) + G1 * Fr(5) == G1 * Fr(12) + G1 * Fr(2)
    assert G1 * Fr(N-1) + G1 == I1


class Fa:
    # A class for elements in polynomial extension fields

    degree = 0
    p = []

    def __init__(self, coeffs):
        assert len(coeffs) == self.degree
        self.coeffs = coeffs

    def __repr__(self):
        return f'Fa({self.coeffs})'

    def __eq__(self, other):
        return self.coeffs == other.coeffs

    def __add__(self, other):
        return self.__class__(polynomial.ext(polynomial.add(self.coeffs, other.coeffs), self.degree))

    def __sub__(self, other):
        return self.__class__(polynomial.ext(polynomial.sub(self.coeffs, other.coeffs), self.degree))

    def __mul__(self, other):
        mulmod = polynomial.rem(polynomial.mul(self.coeffs, other.coeffs), self.p)
        return self.__class__(polynomial.ext(mulmod, self.degree))

    def __truediv__(self, other):
        return self * self.__class__(polynomial.ext(polynomial.inv(other.coeffs, self.p), self.degree))

    def __pow__(self, data):
        result = self.one()
        mulend = self
        while data:
            b = data & 1
            if b == 1:
                result *= mulend
            mulend *= mulend
            data = data >> 1
        return result

    def __pos__(self):
        return self

    def __neg__(self):
        return self.__class__([-c for c in self.coeffs])

    @classmethod
    def nil(cls):
        return cls([Fq(0) for _ in range(cls.degree)])

    @classmethod
    def one(cls):
        return cls([Fq(1)] + [Fq(0) for _ in range(cls.degree - 1)])


class F2(Fa):
    degree = 2
    p = [Fq(e) for e in [1, 0, 1]]  # i² + 1 = 0


if __name__ == '__main__':
    a = F2([Fq(1), Fq(0)])
    b = F2([Fq(1), Fq(2)])
    assert a + b == F2([Fq(2), Fq(2)])
    assert b / b == F2([Fq(1), Fq(0)])
    assert a / b + a / b == (a + a) / b
    assert a * b + a * b == (a + a) * b
    assert a ** (P ** 2 - 1) == a


class Ft(Fa):
    degree = 12
    p = [Fq(e) for e in [82, 0, 0, 0, 0, 0, -18, 0, 0, 0, 0, 0, 1]]  # w¹² - 18w⁶ + 82 = 0


class P2(Pa):
    a = F2([Fq(0), Fq(0)])
    b = F2([Fq(3), Fq(0)]) / F2([Fq(9), Fq(1)])
    i = [
        F2([Fq(0), Fq(0)]),
        F2([Fq(0), Fq(0)])
    ]


class Pt(Pa):
    a = Ft([Fq(0) for _ in range(12)])
    b = Ft([Fq(e) for e in [3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    i = [
        Ft([Fq(0) for _ in range(12)]),
        Ft([Fq(0) for _ in range(12)])
    ]


G2 = P2(
    F2([Fq(0x1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed),
        Fq(0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2)]),
    F2([Fq(0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa),
        Fq(0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b)])
)
I2 = P2(F2.nil(), F2.nil())

if __name__ == '__main__':
    assert G2 * Fr(2) + G2 + G2 == G2 * Fr(4)
    assert G2 + G2 != G2
    assert G2 * Fr(9) + G2 * Fr(5) == G2 * Fr(12) + G2 * Fr(2)
    assert G2 * Fr(N-1) + G2 == I2


def pairing_twist(p):
    if p.x == P2.i[0] and p.y == P2.i[1]:
        return Pt(Pt.i[0], Pt.i[1])
    w = Ft([Fq(e) for e in [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    xcoeffs = [p.x.coeffs[0] - p.x.coeffs[1] * Fq(9), p.x.coeffs[1]]
    ycoeffs = [p.y.coeffs[0] - p.y.coeffs[1] * Fq(9), p.y.coeffs[1]]
    nx = Ft([xcoeffs[0], Fq(0), Fq(0), Fq(0), Fq(0), Fq(0), xcoeffs[1], Fq(0), Fq(0), Fq(0), Fq(0), Fq(0)])
    ny = Ft([ycoeffs[0], Fq(0), Fq(0), Fq(0), Fq(0), Fq(0), ycoeffs[1], Fq(0), Fq(0), Fq(0), Fq(0), Fq(0)])
    return Pt(nx * w ** 2, ny * w ** 3)


Gt = pairing_twist(G2)
It = Pt(Ft.nil(), Ft.nil())

if __name__ == '__main__':
    assert Gt * Fr(2) + Gt + Gt == Gt * Fr(4)
    assert Gt + Gt != Gt
    assert Gt * Fr(9) + Gt * Fr(5) == Gt * Fr(12) + Gt * Fr(2)
    assert Gt * Fr(N-1) + Gt == It


def pairing_line_function(p, q, r):
    # Create a function representing the line between p and q, and evaluate it at r.
    # It can be considered as a distance metric between p + q and the second stationary point r.
    #
    # See https://crypto.stanford.edu/pbc/notes/ep/miller.html
    x1, y1 = p.x, p.y
    x2, y2 = q.x, q.y
    x3, y3 = r.x, r.y
    if x1 != x2:
        m = (y2 - y1) / (x2 - x1)
        return m * (x3 - x1) - (y3 - y1)
    if y1 == y2:
        # Simplify (3x² + a) / 2y => 3x² / 2y
        m = (x1 * x1 + x1 * x1 + x1 * x1) / (y1 + y1)
        return m * (x3 - x1) - (y3 - y1)
    return x3 - x1


if __name__ == '__main__':
    x1, x2, x3 = G1, G1 * Fr(2), G1 * Fr(3)
    y1, y2, y3 = -x1, -x2, -x3
    assert pairing_line_function(x1, x2, x1) == Fq(0)
    assert pairing_line_function(x1, x2, x2) == Fq(0)
    assert pairing_line_function(x1, x2, x3) != Fq(0)
    assert pairing_line_function(x1, x2, y3) == Fq(0)
    assert pairing_line_function(x1, y1, x1) == Fq(0)
    assert pairing_line_function(x1, y1, y1) == Fq(0)
    assert pairing_line_function(x1, y1, x2) != Fq(0)
    assert pairing_line_function(x1, x1, x1) == Fq(0)
    assert pairing_line_function(x1, x1, x2) != Fq(0)
    assert pairing_line_function(x1, x1, y2) == Fq(0)


def pairing_miller_loop(q, p):
    ate_loop_count = 29793968203157093288
    ate_loop_count_log = 63
    if (q.x == Pt.i[0] and q.y == Pt.i[1]) or (p.x == Pt.i[0] and p.y == Pt.i[1]):
        return Ft.one()
    r = q
    f = Ft.one()
    for i in range(ate_loop_count_log, -1, -1):
        f = f * f * pairing_line_function(r, r, p)
        r = r + r
        if ate_loop_count & (2**i):
            f = f * pairing_line_function(r, q, p)
            r = r + q
    a = Pt(q.x ** P, +q.y ** P)
    b = Pt(a.x ** P, -a.y ** P)
    f = f * pairing_line_function(r, a, p)
    r = r + a
    f = f * pairing_line_function(r, b, p)
    return f ** ((P ** 12 - 1) // N)


def pairing(q, p):
    # Pairing computation
    r = Pt(Ft.nil(), Ft.nil())
    r.x.coeffs[0] = p.x
    r.y.coeffs[0] = p.y
    return pairing_miller_loop(pairing_twist(q), r)


if __name__ == '__main__':
    a = pairing(G2, +G1)
    b = pairing(G2, -G1)
    assert a * b == Ft.one()
    c = pairing(-G2, G1)
    assert a * c == Ft.one()
    assert b == c
    assert a ** N == Ft.one()
    d = pairing(G2, G1 * Fr(2))
    assert a * a == d
    assert a != d and a != c and d != c
    e = pairing(G2 * Fr(2), G1)
    assert a * a == e
    f = pairing(G2 * Fr(27), G1 * Fr(37))
    g = pairing(G2, G1 * Fr(999))
    assert f == g
