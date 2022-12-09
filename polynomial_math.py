def polyclr(c1):
    for i in range(len(c1) - 1, -1, -1):
        if c1[i] != 0:
            break
    return c1[:i+1]


def polydeg(c1):
    d = len(c1) - 1
    while c1[d] == 0 and d:
        d -= 1
    return d


def polyadd(c1, c2):
    p = [0 for _ in range(max(len(c1), len(c2)))]
    for i, e in enumerate(c1):
        p[i] += e
    for i, e in enumerate(c2):
        p[i] += e
    return polyclr(p)


def polysub(c1, c2):
    p = [0 for _ in range(max(len(c1), len(c2)))]
    for i, e in enumerate(c1):
        p[i] += e
    for i, e in enumerate(c2):
        p[i] -= e
    return polyclr(p)


def polymul(c1, c2):
    p = [0 for _ in range(len(c1) + len(c2) - 1)]
    for i in range(len(c1)):
        for j in range(len(c2)):
            p[i+j] += c1[i] * c2[j]
    return polyclr(p)


def polydivmod(c1, c2):
    # Algorithm: https://en.wikipedia.org/wiki/Polynomial_long_division
    # The code implementation is inspired by numpy.polynomial.polynomial.polydiv
    lc1 = len(c1)
    lc2 = len(c2)
    if c2[-1] == 0:
        raise ZeroDivisionError()
    if lc1 < lc2:
        return [0], c1
    if lc2 == 1:
        return [e / c2[0] for e in c1], [0]
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
    newt, t = [1], [0]
    newr, r = c1, c2
    while polydeg(newr):
        quotient = polydiv(r, newr)
        r, newr = newr, polysub(r, polymul(newr, quotient))
        t, newt = newt, polysub(t, polymul(newt, quotient))
    return polyclr([e/newr[0] for e in newt[:polydeg(c2)]])


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
assert polymod(polymul(polyinv(px, qx), px), qx)[0] == 1
