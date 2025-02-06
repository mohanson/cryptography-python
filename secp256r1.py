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

    def sqrt(self):
        # https://www.staff.uni-mainz.de/pommeren/Cryptology/Asymmetric/5_NTh/SqRprim.pdf, 5.3
        # Note source number should be quadratic residued.
        if (self.p - 3) % 4 == 0:
            m = (self.p - 3) // 4
            return self ** (m + 1)
        raise Exception('unreachable')

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
P = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
# The order n of G.
N = 0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551


class Fq(Fp):

    p = P

    def __repr__(self):
        return f'Fq(0x{self.x:064x})'


class Fr(Fp):

    p = N

    def __repr__(self):
        return f'Fr(0x{self.x:064x})'


A = Fq(0xffffffff00000001000000000000000000000000fffffffffffffffffffffffc)
B = Fq(0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b)


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
        x1, x2 = self.x, data.x
        y1, y2 = self.y, data.y
        if x1 == Fq(0) and y1 == Fq(0):
            return data
        if x2 == Fq(0) and y2 == Fq(0):
            return self
        if x1 == x2 and y1 == +y2:
            sk = (x1 * x1 + x1 * x1 + x1 * x1 + A) / (y1 + y1)
            x3 = sk * sk - x1 - x2
            y3 = sk * (x1 - x3) - y1
            return Pt(x3, y3)
        if x1 == x2 and y1 == -y2:
            return I
        sk = (y2 - y1) / (x2 - x1)
        x3 = sk * sk - x1 - x2
        y3 = sk * (x1 - x3) - y1
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

    def __truediv__(self, k):
        return self.__mul__(k ** -1)

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
    Fq(0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296),
    Fq(0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5),
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

    p = G * Fr(0x5f6717883bef25f45a129c11fcac1567d74bda5a9ad4cbffc8203c0da2a1473c)
    assert p.x.x == 0x63983e4c8002f443ccb58f7cd8232b75af26c432e30cb7584bed0dbc35bcf86a
    assert p.y.x == 0xcdc066239b4c9a967ffd2429d6ffe57850122163413348ba520726e5b08a9d79

    x = Fq(0x660fe3dd941bc58104fff3b424d82cd69658191f91166af80528e65d07cec0c0)
    assert x.sqrt() * x.sqrt() == x
