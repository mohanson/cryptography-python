P = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f
N = 0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141
G_X = 0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798
G_Y = 0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8
A = 0
B = 7
assert(A < P)
assert(B < P)
assert((4 * A**3 + 27 * B**2) % P != 0)


class Fg:
    # Galois field. In mathematics, a finite field or Galois field is a field that contains a finite number of elements.
    # As with any field, a finite field is a set on which the operations of multiplication, addition, subtraction and
    # division are defined and satisfy certain basic rules.

    p = 0

    def __init__(self, x):
        self.x = x % self.p

    def __repr__(self):
        return f'Fp(0x{self.x:064x})'

    def __eq__(self, other):
        return self.x == other.x

    def __add__(self, other):
        return self.__class__((self.x + other.x) % self.p)

    def __sub__(self, other):
        return self.__class__((self.x - other.x) % self.p)

    def __mul__(self, other):
        return self.__class__((self.x * other.x) % self.p)

    def __truediv__(self, other):
        return self * other ** -1

    def __pow__(self, other):
        return self.__class__(pow(self.x, other, self.p))

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
    def __init__(self, x, y):
        if x != Fp(0) or y != Fp(0):
            assert(y ** 2 == x ** 3 + Fp(A) * x + Fp(B))
        self.x = x
        self.y = y

    def __repr__(self):
        return f'Ec({self.x}, {self.y})'

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def __add__(self, other):
        if self == I:
            return other
        if other == I:
            return self
        if self.x == other.x and self.y == -other.y:
            return I
        x1, x2 = self.x, other.x
        y1, y2 = self.y, other.y
        if self.y == other.y:
            s = (Fp(3) * x1 * x1 + Fp(A)) / (Fp(2) * y1)
        else:
            s = (y2 - y1) / (x2 - x1)
        x3 = s * s - x1 - x2
        y3 = s * (x1 - x3) - y1
        return Ec(x3, y3)

    def __mul__(self, k):
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


# Identity element
I = Ec(Fp(0x0), Fp(0x0))
# Generator point
G = Ec(Fp(G_X), Fp(G_Y))
