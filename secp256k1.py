class Fp:
    # Galois field. In mathematics, a finite field or Galois field is a field that contains a finite number of elements.
    # As with any field, a finite field is a set on which the operations of multiplication, addition, subtraction and
    # division are defined and satisfy certain basic rules.
    #
    # https://www.cs.miami.edu/home/burt/learning/Csc609.142/ecdsa-cert.pdf
    # Don Johnson, Alfred Menezes and Scott Vanstone, The Elliptic Curve Digital Signature Algorithm (ECDSA)
    # 3.1 The Finite Field Fp

    def __init__(self, p, x):
        self.p = p
        self.x = x % self.p

    def __repr__(self):
        return f'Fp(0x{self.x:064x})'

    def __eq__(self, data):
        assert self.p == data.p
        return self.x == data.x

    def __add__(self, data):
        assert self.p == data.p
        return Fp(self.p, (self.x + data.x) % self.p)

    def __sub__(self, data):
        assert self.p == data.p
        return Fp(self.p, (self.x - data.x) % self.p)

    def __mul__(self, data):
        assert self.p == data.p
        return Fp(self.p, (self.x * data.x) % self.p)

    def __truediv__(self, data):
        assert self.p == data.p
        return self * data ** -1

    def __pow__(self, data):
        return Fp(self.p, pow(self.x, data, self.p))

    def __neg__(self):
        return Fp(self.p, self.p - self.x)


if __name__ == '__main__':
    assert Fp(23, 12) + Fp(23, 20) == Fp(23, 9)
    assert Fp(23, 8) * Fp(23, 9) == Fp(23, 3)
    assert Fp(23, 8) ** -1 == Fp(23, 3)

P = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f
N = 0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141


class Fq(Fp):

    def __init__(self, x):
        super(Fq, self).__init__(P, x)

    def __repr__(self):
        return f'Fq(0x{self.x:064x})'


class Fr(Fp):

    def __init__(self, x):
        super(Fr, self).__init__(N, x)

    def __repr__(self):
        return f'Fr(0x{self.x:064x})'


G_X = 0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798
G_Y = 0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8
A = 0
B = 7
assert A < P
assert B < P
assert (4 * A**3 + 27 * B**2) % P != 0


class Point:
    a = Fq(A)
    b = Fq(B)
    inf_x = Fq(0)
    inf_y = Fq(0)

    def __init__(self, x, y):
        if x != self.inf_x or y != self.inf_y:
            assert y ** 2 == x ** 3 + self.a * x + self.b
        self.x = x
        self.y = y

    def __repr__(self):
        return f'Point({self.x}, {self.y})'

    def __eq__(self, data):
        return self.x == data.x and self.y == data.y

    def __add__(self, data):
        # https://www.cs.miami.edu/home/burt/learning/Csc609.142/ecdsa-cert.pdf
        # Don Johnson, Alfred Menezes and Scott Vanstone, The Elliptic Curve Digital Signature Algorithm (ECDSA)
        # 4.1 Elliptic Curves Over Fp
        if self.x == self.inf_x and self.y == self.inf_y:
            return data
        if data.x == self.inf_x and data.y == self.inf_y:
            return self
        if self.x == data.x and self.y == -data.y:
            return self.__class__(self.inf_x, self.inf_y)
        x1, x2 = self.x, data.x
        y1, y2 = self.y, data.y
        if self.y == data.y:
            s = (x1 * x1 + x1 * x1 + x1 * x1 + self.a) / (y1 + y1)
        else:
            s = (y2 - y1) / (x2 - x1)
        x3 = s * s - x1 - x2
        y3 = s * (x1 - x3) - y1
        return self.__class__(x3, y3)

    def __sub__(self, data):
        return self + data.__neg__()

    def __mul__(self, k):
        # Point multiplication: Double-and-add
        # https://en.wikipedia.org/wiki/Elliptic_curve_point_multiplication
        n = k.x
        result = self.__class__(self.inf_x, self.inf_y)
        addend = self
        while n:
            b = n & 1
            if b == 1:
                result += addend
            addend = addend + addend
            n = n >> 1
        return result

    def __neg__(self):
        return self.__class__(self.x, -self.y)


# Identity element
I = Point(Fq(0x0), Fq(0x0))
# Generator point
G = Point(Fq(G_X), Fq(G_Y))
