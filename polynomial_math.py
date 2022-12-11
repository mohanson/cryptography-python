class PolynomialCalculator:
    nil = 0
    one = 1

    @classmethod
    def clr(cls, c1):
        for i in range(len(c1) - 1, -1, -1):
            if c1[i] != cls.nil:
                break
        return c1[:i+1]

    @classmethod
    def deg(cls, c1):
        d = len(c1) - 1
        while c1[d] == cls.nil and d:
            d -= 1
        return d

    @classmethod
    def ext(cls, c1, sz):
        p = [cls.nil for _ in range(sz)]
        for i, e in enumerate(c1):
            p[i] = e
        return p

    @classmethod
    def add(cls, c1, c2):
        p = [cls.nil for _ in range(max(len(c1), len(c2)))]
        for i, e in enumerate(c1):
            p[i] += e
        for i, e in enumerate(c2):
            p[i] += e
        return cls.clr(p)

    @classmethod
    def sub(cls, c1, c2):
        p = [cls.nil for _ in range(max(len(c1), len(c2)))]
        for i, e in enumerate(c1):
            p[i] += e
        for i, e in enumerate(c2):
            p[i] -= e
        return cls.clr(p)

    @classmethod
    def mul(cls, c1, c2):
        p = [cls.nil for _ in range(len(c1) + len(c2) - 1)]
        for i in range(len(c1)):
            for j in range(len(c2)):
                p[i+j] += c1[i] * c2[j]
        return cls.clr(p)

    @classmethod
    def divrem(cls, c1, c2):
        # Algorithm: https://en.wikipedia.org/wiki/Polynomial_long_division
        # The code implementation is inspired by numpy.polynomial.polynomial.polydiv
        lc1 = len(c1)
        lc2 = len(c2)
        if c2[-1] == cls.nil:
            raise ZeroDivisionError()
        if lc1 < lc2:
            return [cls.nil], c1
        if lc2 == 1:
            return [e / c2[0] for e in c1], [cls.nil]
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
        return [e/scl for e in nc1[j+1:]], cls.clr(nc1[:j+1])

    @classmethod
    def div(cls, c1, c2):
        return cls.divrem(c1, c2)[0]

    @classmethod
    def rem(cls, c1, c2):
        return cls.divrem(c1, c2)[1]

    @classmethod
    def inv(cls, c1, c2):
        # Algorithm: https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
        newt, t = [cls.one], [cls.nil]
        newr, r = c1, c2
        while cls.deg(newr):
            quotient = cls.div(r, newr)
            r, newr = newr, cls.sub(r, cls.mul(newr, quotient))
            t, newt = newt, cls.sub(t, cls.mul(newt, quotient))
        return cls.clr([e/newr[0] for e in newt[:cls.deg(c2)]])


if __name__ == '__main__':
    px = [4, -2, 5]
    qx = [2, -5, 2]
    assert PolynomialCalculator.add(px, qx) == [6, -7, 7]
    assert PolynomialCalculator.sub(px, qx) == [2, 3, 3]
    assert PolynomialCalculator.mul(px, qx) == [8, -24, 28, -29, 10]
    assert PolynomialCalculator.div(px, qx) == [2.5]
    assert PolynomialCalculator.rem(px, qx) == [-1, 10.5]
    assert PolynomialCalculator.sub(px, px) == [0]
    # https://en.wikipedia.org/wiki/Polynomial_long_division#Example
    assert PolynomialCalculator.div([-4, 0, -2, 1], [-3, 1]) == [3, 1, 1]
    assert PolynomialCalculator.rem([-4, 0, -2, 1], [-3, 1]) == [5]
    assert PolynomialCalculator.rem(PolynomialCalculator.mul(PolynomialCalculator.inv(px, qx), px), qx)[0] == 1
