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
