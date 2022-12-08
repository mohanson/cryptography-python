def polyadd(c1, c2):
    return [c1[i] + c2[i] for i in range(len(c1))]


def polysub(c1, c2):
    return [c1[i] - c2[i] for i in range(len(c1))]


def polymul(c1, c2):
    p = [0 for _ in range(len(c1) + len(c2) - 1)]
    for i in range(len(c1)):
        for j in range(len(c2)):
            p[i+j] += c1[i] * c2[j]
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


px = [4, -2, 5]
qx = [2, -5, 2]
assert polyadd(px, qx) == [6, -7, 7]
assert polysub(px, qx) == [2, 3, 3]
assert polymul(px, qx) == [8, -24, 28, -29, 10]
assert polydiv(px, qx) == [2.5]
assert polymod(px, qx) == [-1, 10.5]
# https://en.wikipedia.org/wiki/Polynomial_long_division#Example
assert polydiv([-4, 0, -2, 1], [-3, 1]) == [3, 1, 1]
assert polymod([-4, 0, -2, 1], [-3, 1]) == [5]
