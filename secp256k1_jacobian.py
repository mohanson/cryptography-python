import secp256k1

# [1] 熠智科技, 椭圆曲线科普, https://download.yeez.tech/doc/ECcurve.pdf, 3.4 Projective Space
# [2] https://en.wikibooks.org/wiki/Cryptography/Prime_Curve/Jacobian_Coordinates


class Pj:

    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return f'Pj({self.x}, {self.y}, {self.z})'

    def __eq__(self, data):
        if self.x * data.z != data.x * self.z:
            return 0
        if self.y * data.z != data.y * self.z:
            return 0
        return 1

    @classmethod
    def encode(cls, pt):
        if pt == secp256k1.I:
            return Pj(secp256k1.Fq(0), secp256k1.Fq(1), secp256k1.Fq(0))
        else:
            return Pj(pt.x, pt.y, secp256k1.Fq(1))

    def decode(self):
        if self.z == secp256k1.Fq(0):
            return secp256k1.I
        else:
            return secp256k1.Pt(self.x / self.z, self.y / self.z)

    def __add__(self, data):
        x1, y1, z1 = self.x, self.y, self.z
        x2, y2, z2 = data.x, data.y, data.z
        if z1 == secp256k1.Fq(0):
            return data
        if z2 == secp256k1.Fq(0):
            return self
        u1 = y2 * z1
        u2 = y1 * z2
        v1 = x2 * z1
        v2 = x1 * z2
        if v1 == v2:
            if u1 != u2:
                return I
            else:
                t = secp256k1.A * z1 * z1 + secp256k1.Fq(3) * x1 * x1
                u = y1 * z1
                v = u * x1 * y1
                w = t * t - secp256k1.Fq(8) * v
                x3 = secp256k1.Fq(2) * u * w
                y3 = t * (secp256k1.Fq(4) * v - w) - secp256k1.Fq(8) * y1 * y1 * u * u
                z3 = secp256k1.Fq(8) * u * u * u
                return Pj(x3, y3, z3)
        else:
            u = u1 - u2
            v = v1 - v2
            w = u * u * z1 * z2 - v * v * v - secp256k1.Fq(2) * v * v * x1 * z2
            x3 = v * w
            y3 = u * (v * v * x1 * z2 - w) - v * v * v * y1 * z2
            z3 = v * v * v * z1 * z2
            return Pj(x3, y3, z3)

    def __sub__(self, data):
        return self + data.__neg__()

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

    def __neg__(self):
        return Pj(self.x, -self.y, self.z)


G = Pj.encode(secp256k1.G)
I = Pj.encode(secp256k1.I)


if __name__ == '__main__':
    p1 = secp256k1.G * secp256k1.Fr(42)
    p2 = secp256k1.G * secp256k1.Fr(24)
    j1 = Pj.encode(p1)
    j2 = Pj.encode(p2)
    assert p1 == j1.decode()
    assert p2 == j2.decode()
    assert j1 == G * secp256k1.Fr(42)
    assert j2 == G * secp256k1.Fr(24)
    assert p1 + p2 == (j1 + j2).decode()
    assert p1 + p1 == (j1 + j1).decode()
    assert j1 - j1 == I
    assert j1 + I == j1
    assert j1 * secp256k1.Fr(42) == Pj.encode(p1 * secp256k1.Fr(42))
