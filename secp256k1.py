P = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f
N = 0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141
G_X = 0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798
G_Y = 0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8
A = 0
B = 7
assert A < P
assert B < P
assert (4 * A**3 + 27 * B**2) % P != 0


class Fg:
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
        return self.x == data.x

    def __add__(self, data):
        return self.__class__((self.x + data.x) % self.p)

    def __sub__(self, data):
        return self.__class__((self.x - data.x) % self.p)

    def __mul__(self, data):
        return self.__class__((self.x * data.x) % self.p)

    def __truediv__(self, data):
        return self * data ** -1

    def __pow__(self, data):
        return self.__class__(pow(self.x, data, self.p))

    def __neg__(self):
        return self.__class__(self.p - self.x)


Fg.p = 23
assert Fg(12) + Fg(20) == Fg(9)
assert Fg(8) * Fg(9) == Fg(3)
assert Fg(8) ** -1 == Fg(3)


class Fp(Fg):
    p = P

    def __repr__(self):
        return f'Fp(0x{self.x:064x})'


class Fr(Fg):
    p = N

    def __repr__(self):
        return f'Fr(0x{self.x:064x})'


class Ec:
    a = Fp(A)
    b = Fp(B)
    inf_x = Fp(0)
    inf_y = Fp(0)

    def __init__(self, x, y):
        if x != self.inf_x or y != self.inf_y:
            assert y ** 2 == x ** 3 + self.a * x + self.b
        self.x = x
        self.y = y

    def __repr__(self):
        return f'Ec({self.x}, {self.y})'

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
            return I
        x1, x2 = self.x, data.x
        y1, y2 = self.y, data.y
        if self.y == data.y:
            s = (x1 * x1 + x1 * x1 + x1 * x1 + self.a) / (y1 + y1)
        else:
            s = (y2 - y1) / (x2 - x1)
        x3 = s * s - x1 - x2
        y3 = s * (x1 - x3) - y1
        return self.__class__(x3, y3)

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

    def __neg__(self):
        return self.__class__(self.x, -self.y)


# Identity element
I = Ec(Fp(0x0), Fp(0x0))
# Generator point
G = Ec(Fp(G_X), Fp(G_Y))
