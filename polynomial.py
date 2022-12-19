def clr(c1):
    for i in range(len(c1) - 1, -1, -1):
        if c1[i] != c1[0].__class__(0):
            break
    return c1[:i+1]


def deg(c1):
    d = len(c1) - 1
    while c1[d] == c1[0].__class__(0) and d:
        d -= 1
    return d


def ext(c1, sz):
    p = [c1[0].__class__(0) for _ in range(sz)]
    for i, e in enumerate(c1):
        p[i] = e
    return p


def add(c1, c2):
    p = [c1[0].__class__(0) for _ in range(max(len(c1), len(c2)))]
    for i, e in enumerate(c1):
        p[i] += e
    for i, e in enumerate(c2):
        p[i] += e
    return clr(p)


def sub(c1, c2):
    p = [c1[0].__class__(0) for _ in range(max(len(c1), len(c2)))]
    for i, e in enumerate(c1):
        p[i] += e
    for i, e in enumerate(c2):
        p[i] -= e
    return clr(p)


def mul(c1, c2):
    p = [c1[0].__class__(0) for _ in range(len(c1) + len(c2) - 1)]
    for i in range(len(c1)):
        for j in range(len(c2)):
            p[i+j] += c1[i] * c2[j]
    return clr(p)


def divrem(c1, c2):
    # Algorithm: https://en.wikipedia.org/wiki/Polynomial_long_division
    # The code implementation is inspired by numpy.polynomial.polynomial.polydiv
    lc1 = len(c1)
    lc2 = len(c2)
    if c2[-1] == c1[0].__class__(0):
        raise ZeroDivisionError()
    if lc1 < lc2:
        return [c1[0].__class__(0)], c1
    if lc2 == 1:
        return [e / c2[0] for e in c1], [c1[0].__class__(0)]
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
    return [e/scl for e in nc1[j+1:]], clr(nc1[:j+1])


def div(c1, c2):
    return divrem(c1, c2)[0]


def rem(c1, c2):
    return divrem(c1, c2)[1]


def inv(c1, c2):
    # Algorithm: https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
    newt, t = [c1[0].__class__(1)], [c1[0].__class__(0)]
    newr, r = c1, c2
    while deg(newr):
        quotient = div(r, clr(newr))
        r, newr = newr, sub(r, mul(newr, quotient))
        t, newt = newt, sub(t, mul(newt, quotient))
    return clr([e/newr[0] for e in newt[:deg(c2)]])


def interp(c1, c2):
    # Lagrange interpolation, copied from scipy.interpolate.lagrange.
    M = len(c1)
    p = [c1[0].__class__(0)]
    for j in range(M):
        pt = [c2[j]]
        for k in range(M):
            if k == j:
                continue
            fac = c1[j]-c1[k]
            pt = mul([c1[0].__class__(1) / fac, -c1[k] / fac], pt)
        p = add(p, pt)
    return p


def vanish(c1):
    # See: https://aszepieniec.github.io/stark-anatomy/basic-tools
    x = [c1[0].__class__(0), c1[0].__class__(1)]
    a = [c1[0].__class__(1)]
    for d in c1:
        a = mul(a, sub(x, [d]))
    return a


if __name__ == '__main__':
    c1 = [4, -2, 5]
    c2 = [2, -5, 2]
    assert add(c1, c2) == [6, -7, 7]
    assert sub(c1, c2) == [2, 3, 3]
    assert mul(c1, c2) == [8, -24, 28, -29, 10]
    assert div(c1, c2) == [2.5]
    assert rem(c1, c2) == [-1, 10.5]
    assert sub(c1, c1) == [0]
    # Copied from https://en.wikipedia.org/wiki/Polynomial_long_division#Example
    assert div([-4, 0, -2, 1], [-3, 1]) == [3, 1, 1]
    assert rem([-4, 0, -2, 1], [-3, 1]) == [5]
    assert rem(mul(inv(c1, c2), c1), c2)[0] == 1
    assert interp([1,  2,  3,  4], [4, 15, 40, 85]) == [0.9999999999999982, 1.0, 1.0000000000000284, 1.0]
    assert vanish([0, 1]) == [0, -1, 1]
    assert vanish([1, 2]) == [2, -3, 1]
