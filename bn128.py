import secp256k1
import polynomial_math
import sys
sys.setrecursionlimit(10000)


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


class Ec(secp256k1.Ec):
    a = Fp(0)
    b = Fp(3)
    inf_x = Fp(0)
    inf_y = Fp(0)


class Ec2(secp256k1.Ec):
    a = Fp2([Fp(0), Fp(0)])
    b = Fp2([Fp(3), Fp(0)]) / Fp2([Fp(9), Fp(1)])
    inf_x = Fp2([Fp(0), Fp(0)])
    inf_y = Fp2([Fp(0), Fp(0)])

    def twist(self):
        if self.x == self.inf_x and self.y == self.inf_y:
            return Ec12(Ec12.inf_x, Ec12.inf_y)
        # "Twist" a point in E(FQ2) into a point in E(FQ12)
        w = Fp12([Fp(e) for e in [0, 1] + [0] * 10])
        # Field isomorphism from Z[p] / x**2 to Z[p] / x**2 - 18*x + 82
        xcoeffs = [self.x.coeffs[0] - self.x.coeffs[1] * Fp(9), self.x.coeffs[1]]
        ycoeffs = [self.y.coeffs[0] - self.y.coeffs[1] * Fp(9), self.y.coeffs[1]]
        # Isomorphism into subfield of Z[p] / w**12 - 18 * w**6 + 82,
        # where w**6 = x
        nx = Fp12([xcoeffs[0], Fp(0), Fp(0), Fp(0), Fp(0), Fp(0), xcoeffs[1], Fp(0), Fp(0), Fp(0), Fp(0), Fp(0)])
        ny = Fp12([ycoeffs[0], Fp(0), Fp(0), Fp(0), Fp(0), Fp(0), ycoeffs[1], Fp(0), Fp(0), Fp(0), Fp(0), Fp(0)])
        # Divide x coord by w**2 and y coord by w**3
        return Ec12(nx * w ** 2, ny * w**3)


class Ec12(secp256k1.Ec):
    a = Fp12([Fp(0) for _ in range(12)])
    b = Fp12([Fp(e) for e in [3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    inf_x = Fp12([Fp(0) for _ in range(12)])
    inf_y = Fp12([Fp(0) for _ in range(12)])


G1 = Ec(Fp(1), Fp(2))
G2 = Ec2(
    Fp2([Fp(0x1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed),
         Fp(0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2)]),
    Fp2([Fp(0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa),
         Fp(0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b)])
)
G12 = G2.twist()

# Create a function representing the line between P1 and P2,
# and evaluate it at T


def linefunc(P1, P2, T):
    # assert P1 and P2 and T # No points-at-infinity allowed, sorry
    x1, y1 = P1.x, P1.y
    x2, y2 = P2.x, P2.y
    xt, yt = T.x, T.y
    if x1 != x2:
        m = (y2 - y1) / (x2 - x1)
        return m * (xt - x1) - (yt - y1)
    elif y1 == y2:
        m = (x1*x1 + x1*x1 + x1*x1) / (y1 + y1)
        return m * (xt - x1) - (yt - y1)
    else:
        return xt - x1


# Check consistency of the "line function"
one, two, three = G1, G1 * Fr(2), G1 * Fr(3)
negone, negtwo, negthree = G1 * Fr(N - 1), G1 * Fr(N - 2), G1 * Fr(N - 3)

assert linefunc(one, two, one) == Fp(0)
assert linefunc(one, two, two) == Fp(0)
assert linefunc(one, two, three) != Fp(0)
assert linefunc(one, two, negthree) == Fp(0)
assert linefunc(one, negone, one) == Fp(0)
assert linefunc(one, negone, negone) == Fp(0)
assert linefunc(one, negone, two) != Fp(0)
assert linefunc(one, one, one) == Fp(0)
assert linefunc(one, one, two) != Fp(0)
assert linefunc(one, one, negtwo) == Fp(0)


def cast_point_to_fq12(pt):
    if pt.x == Fp(0) and pt.y == Fp(0):
        return Ec12(Fp12([Fp(0) for _ in range(12)]), Fp12([Fp(0) for _ in range(12)]))
    x, y = pt.x, pt.y
    return Ec12(Fp12([x] + [Fp(0)] * 11), Fp12([y] + [Fp(0)] * 11))


ate_loop_count = 29793968203157093288
log_ate_loop_count = 63


def miller_loop(q, p):
    # Main miller loop
    if (q.x == Ec12.inf_x and q.y == Ec12.inf_y) or (p.x == Ec12.inf_x and p.y == Ec12.inf_y):
        return Fp12([Fp(e) for e in [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    R = q
    f = Fp12([Fp(e) for e in [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])  # FQ12.one()
    for i in range(log_ate_loop_count, -1, -1):
        f = f * f * linefunc(R, R, p)
        R = R + R
        if ate_loop_count & (2**i):
            f = f * linefunc(R, q, p)
            R = R + q
    # assert R == multiply(Q, ate_loop_count)
    Q1 = Ec12(q.x ** P, q.y ** P)
    # assert is_on_curve(Q1, b12)
    nQ2 = Ec12(Q1.x ** P, -Q1.y ** P)
    # assert is_on_curve(nQ2, b12)
    f = f * linefunc(R, Q1, p)
    R = R + Q1
    f = f * linefunc(R, nQ2, p)
    # R = add(R, nQ2) This line is in many specifications but it technically does nothing
    return f ** ((P ** 12 - 1) // N)


def pairing(Q, P):
    # Pairing computation
    return miller_loop(Q.twist(), cast_point_to_fq12(P))


p1 = pairing(G2, G1)
pn1 = pairing(G2, -G1)
assert p1 * pn1 == Fp12([Fp(e) for e in [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
np1 = pairing(-G2, G1)
assert p1 * np1 == Fp12([Fp(e) for e in [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
assert pn1 == np1
assert p1 ** N == Fp12([Fp(e) for e in [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
p2 = pairing(G2, G1 * Fr(2))
assert p1 * p1 == p2
assert p1 != p2 and p1 != np1 and p2 != np1
po2 = pairing(G2 * Fr(2), G1)
assert p1 * p1 == po2
p3 = pairing(G2 * Fr(27), G1 * Fr(37))
po3 = pairing(G2, G1 * Fr(999))
assert p3 == po3
