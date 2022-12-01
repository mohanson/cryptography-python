P = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f
N = 0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141
G_X = 0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798
G_Y = 0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8
A = 0
B = 7
assert(A < P)
assert(B < P)
assert((4 * A**3 + 27 * B**2) % P != 0)


class Fp:
    def __init__(self, x):
        self.x = x % P

    def __repr__(self):
        return f'Fp(0x{self.x:064x})'

    def __eq__(self, other):
        return self.x == other.x

    def __add__(self, other):
        return Fp((self.x + other.x) % P)

    def __sub__(self, other):
        return Fp((self.x - other.x) % P)

    def __mul__(self, other):
        return Fp((self.x * other.x) % P)

    def __truediv__(self, other):
        return self * other ** -1

    def __pow__(self, other):
        return Fp(pow(self.x, other, P))

    def __neg__(self):
        return Fp(P - self.x)


class Fr:
    def __init__(self, x):
        self.x = x % N

    def __repr__(self):
        return f'Fr(0x{self.x:064x})'

    def __eq__(self, other):
        return self.x == other.x

    def __add__(self, other):
        return Fr((self.x + other.x) % N)

    def __sub__(self, other):
        return Fr((self.x - other.x) % N)

    def __mul__(self, other):
        return Fr((self.x * other.x) % N)

    def __truediv__(self, other):
        return self * other ** -1

    def __pow__(self, other):
        return Fr(pow(self.x, other, N))

    def __neg__(self):
        return Fr(N - self.x)


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
