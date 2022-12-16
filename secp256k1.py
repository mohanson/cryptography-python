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
P = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f
# The order n of G.
N = 0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141


class Fq(Fp):

    p = P

    def __repr__(self):
        return f'Fq(0x{self.x:064x})'


class Fr(Fp):

    p = N

    def __repr__(self):
        return f'Fr(0x{self.x:064x})'


A = Fq(0)
B = Fq(7)


class Pt:

    def __init__(self, x, y):
        if x != Fq(0) or y != Fq(0):
            assert y ** 2 == x ** 3 + A * x + B
        self.x = x
        self.y = y

    def __repr__(self):
        return f'Pt({self.x}, {self.y})'

    def __eq__(self, data):
        return self.x == data.x and self.y == data.y

    def __add__(self, data):
        # https://www.cs.miami.edu/home/burt/learning/Csc609.142/ecdsa-cert.pdf
        # Don Johnson, Alfred Menezes and Scott Vanstone, The Elliptic Curve Digital Signature Algorithm (ECDSA)
        # 4.1 Elliptic Curves Over Fp
        if self.x == Fq(0) and self.y == Fq(0):
            return data
        if data.x == Fq(0) and data.y == Fq(0):
            return self
        if self.x == data.x and self.y == -data.y:
            return I
        x1, x2 = self.x, data.x
        y1, y2 = self.y, data.y
        if self.y == data.y:
            s = (x1 * x1 + x1 * x1 + x1 * x1 + A) / (y1 + y1)
        else:
            s = (y2 - y1) / (x2 - x1)
        x3 = s * s - x1 - x2
        y3 = s * (x1 - x3) - y1
        return Pt(x3, y3)

    def __sub__(self, data):
        return self + data.__neg__()

    def __mul__(self, k):
        # Point multiplication: Double-and-add
        # https://en.wikipedia.org/wiki/Elliptic_curve_point_multiplication
        n = k.x
        result = I
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
        return Pt(self.x, -self.y)


# Identity element
I = Pt(
    Fq(0),
    Fq(0),
)
# Generator point
G = Pt(
    Fq(0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798),
    Fq(0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8),
)

if __name__ == '__main__':
    p = G * Fr(42)
    q = G * Fr(24)
    r = Pt(p.x, -p.y)
    assert p + q == G * Fr(66)
    assert p + p == G * Fr(84)
    assert p - q == G * Fr(18)
    assert r == -p
    assert p + r == I
    assert p + I == p
    assert p * Fr(42) == G * Fr(1764)
