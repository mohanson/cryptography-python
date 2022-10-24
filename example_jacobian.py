import secp256k1


class EcJacobian:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    @classmethod
    def encode(cls, ec):
        if ec == secp256k1.I:
            return EcJacobian(secp256k1.Fp(0), secp256k1.Fp(1), secp256k1.Fp(0))
        else:
            return EcJacobian(ec.x, ec.y, secp256k1.Fp(1))

    def decode(self):
        if self.z == secp256k1.Fp(0):
            return secp256k1.I
        else:
            return secp256k1.Ec(self.x / self.z, self.y / self.z)

    def double(self):
        x, y, z = self.x, self.y, self.z
        a = secp256k1.Fp(3) * x * x
        b = y * z
        c = x * y * b
        d = a * a - secp256k1.Fp(8) * c
        e = b * b
        f = secp256k1.Fp(2) * d * b
        g = a * (secp256k1.Fp(4) * c - d) - secp256k1.Fp(8) * y * y * e
        h = secp256k1.Fp(8) * b * e
        return EcJacobian(f, g, h)

    def __add__(self, other):
        x1, y1, z1 = self.x, self.y, self.z
        x2, y2, z2 = other.x, other.y, other.z
        if z1 == secp256k1.Fp(0):
            return other
        if z2 == secp256k1.Fp(0):
            return self
        u1 = y2 * z1
        u2 = y1 * z2
        v1 = x2 * z1
        v2 = x1 * z2
        if v1 == v2:
            if u1 != u2:
                return EcJacobian.encode(secp256k1.I)
            else:
                return self.double()
        u = u1 - u2
        v = v1 - v2
        w = z1 * z2
        a = v * v
        b = a * v2
        c = v * a
        d = u * u * w - c - b * secp256k1.Fp(2)
        x3 = v * d
        y3 = u * (b - d) - c * u2
        z3 = c * w
        return EcJacobian(x3, y3, z3)


p = secp256k1.G * secp256k1.Fr(42)
q = secp256k1.G * secp256k1.Fr(24)
assert p + q == secp256k1.G * secp256k1.Fr(66)
assert (EcJacobian.encode(p) + EcJacobian.encode(q)).decode() == secp256k1.G * secp256k1.Fr(66)

assert p + p == secp256k1.G * secp256k1.Fr(84)
assert (EcJacobian.encode(p) + EcJacobian.encode(p)).decode() == secp256k1.G * secp256k1.Fr(84)

q = secp256k1.Ec(p.x, -p.y)
assert p + q == secp256k1.I
assert (EcJacobian.encode(p) + EcJacobian.encode(q)).decode() == secp256k1.I
