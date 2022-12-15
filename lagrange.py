import secp256k1
import polynomial_math

poly = polynomial_math.PolynomialCalculator


def lagrange(x, w):
    # Copied from scipy.interpolate.lagrange
    M = len(x)
    p = [0.0]
    for j in range(M):
        pt = [w[j]]
        for k in range(M):
            if k == j:
                continue
            fac = x[j]-x[k]
            pt = poly.mul([1.0 / fac, -x[k] / fac], pt)
        p = poly.add(p, pt)
    return p


x = [1,  2,  3,  4]
y = [4, 15, 40, 85]
ret = lagrange(x, y)
assert ret == [0.9999999999999982, 1.0, 1.0000000000000284, 1.0]


class Fq(secp256k1.Fp):

    def __init__(self, x):
        super(Fq, self).__init__(13, x)


class Foly(poly):
    nil = Fq(0)
    one = Fq(1)


def lagrange_foly(x, w):
    # Copied from scipy.interpolate.lagrange
    M = len(x)
    p = [Foly.nil]
    for j in range(M):
        pt = [w[j]]
        for k in range(M):
            if k == j:
                continue
            fac = x[j]-x[k]
            pt = Foly.mul([Foly.one / fac, -x[k] / fac], pt)
        p = Foly.add(p, pt)
    return p


assert lagrange_foly([Fq(1), Fq(4)], [Fq(6), Fq(2)]) == [Fq(3), Fq(3)]  # 3x + 3
