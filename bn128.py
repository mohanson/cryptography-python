import polynomial_math
import secp256k1

P = 0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47
N = 0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001

assert pow(2, N, N) == 2
assert (P ** 12 - 1) % N == 0


class Fp(secp256k1.Fg):
    p = P

    def __repr__(self):
        return f'Fp(0x{self.x:064x})'


class Fr(secp256k1.Fg):
    p = N

    def __repr__(self):
        return f'Fr(0x{self.x:064x})'


class Pc(polynomial_math.PolynomialCalculator):
    nil = Fp(0)
    one = Fp(1)


class Fpx:
    # A class for elements in polynomial extension fields

    degree = 0
    p = []
    nil = Fp(0)
    one = Fp(1)

    def __init__(self, coeffs):
        assert len(coeffs) == self.degree
        self.coeffs = coeffs

    def __repr__(self):
        return f'Fpx({self.coeffs})'

    def __eq__(self, other):
        return self.coeffs == other.coeffs

    def __add__(self, other):
        return self.__class__(Pc.ext(Pc.add(self.coeffs, other.coeffs), self.degree))

    def __sub__(self, other):
        return self.__class__(Pc.ext(Pc.sub(self.coeffs, other.coeffs), self.degree))

    def __mul__(self, other):
        return self.__class__(Pc.ext(Pc.rem(Pc.mul(self.coeffs, other.coeffs), self.p), self.degree))

    def __truediv__(self, other):
        return self * self.__class__(Pc.ext(Pc.inv(other.coeffs, self.p), self.degree))

    def __pow__(self, other):
        if other == 0:
            return self.__class__([self.one] + [self.nil for _ in range(self.degree - 1)])
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


# Curve is y**2 = x**3 + 3
B = Fp(3)
# Twisted curve over FQ**2
B2 = Fp2([Fp(3), Fp(0)]) / Fp2([Fp(9), Fp(1)])
# Extension curve over FQ**12; same b value as over FQ
B12 = Fp12([Fp(e) for e in [3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])


# Generator for curve over FQ
G1 = (Fp(1), Fp(2))
# Generator for twisted curve over FQ2
G2 = (Fp2([Fp(0x1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed), Fp(0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2)]),
      Fp2([Fp(0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa), Fp(0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b)]))


def is_inf(pt):
    # Check if a point is the point at infinity
    return pt is None


def is_on_curve(pt, b):
    # Check that a point is on the curve defined by y**2 == x**3 + b
    if is_inf(pt):
        return True
    x, y = pt
    return y**2 - x**3 == b


assert is_on_curve(G1, B)
assert is_on_curve(G2, B2)


def double(pt):
    # Elliptic curve doubling
    x, y = pt
    l = 3 * x**2 / (2 * y)
    newx = l**2 - 2 * x
    newy = -l * newx + l * x - y
    return newx, newy


def add(p1, p2):
    # Elliptic curve addition
    if p1 is None or p2 is None:
        return p1 if p2 is None else p2
    x1, y1 = p1
    x2, y2 = p2
    if x2 == x1 and y2 == y1:
        return double(p1)
    elif x2 == x1:
        return None
    else:
        l = (y2 - y1) / (x2 - x1)
    newx = l**2 - x1 - x2
    newy = -l * newx + l * x1 - y1
    assert newy == (-l * newx + l * x2 - y2)
    return (newx, newy)


def multiply(pt, n):
    # Elliptic curve point multiplication
    if n == 0:
        return None
    elif n == 1:
        return pt
    elif not n % 2:
        return multiply(double(pt), n // 2)
    else:
        return add(multiply(double(pt), int(n // 2)), pt)


def eq(p1, p2):
    return p1 == p2


# "Twist" a point in E(FQ2) into a point in E(FQ12)
w = Fp12([Fp(e) for e in [0, 1] + [0] * 10])


def neg(pt):
    # Convert P => -P
    if pt is None:
        return None
    x, y = pt
    return (x, -y)


def twist(pt):
    if pt is None:
        return None
    _x, _y = pt
    # Field isomorphism from Z[p] / x**2 to Z[p] / x**2 - 18*x + 82
    xcoeffs = [_x.coeffs[0] - _x.coeffs[1] * Fp(9), _x.coeffs[1]]
    ycoeffs = [_y.coeffs[0] - _y.coeffs[1] * Fp(9), _y.coeffs[1]]
    # Isomorphism into subfield of Z[p] / w**12 - 18 * w**6 + 82,
    # where w**6 = x
    nx = Fp12([xcoeffs[0], Fp(0), Fp(0), Fp(0), Fp(0), Fp(0), xcoeffs[1], Fp(0), Fp(0), Fp(0), Fp(0), Fp(0)])
    ny = Fp12([ycoeffs[0], Fp(0), Fp(0), Fp(0), Fp(0), Fp(0), ycoeffs[1], Fp(0), Fp(0), Fp(0), Fp(0), Fp(0)])
    # Divide x coord by w**2 and y coord by w**3
    return (nx * w ** 2, ny * w**3)


G12 = twist(G2)
# Check that the twist creates a point that is on the curve
assert is_on_curve(G12, B12)
