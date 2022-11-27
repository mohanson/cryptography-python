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
                t = secp256k1.Fp(secp256k1.A) * z1 * z1 + secp256k1.Fp(3) * x1 * x1
                u = y1 * z1
                v = u * x1 * y1
                w = t * t - secp256k1.Fp(8) * v
                x3 = secp256k1.Fp(2) * u * w
                y3 = t * (secp256k1.Fp(4) * v - w) - secp256k1.Fp(8) * y1 * y1 * u * u
                z3 = secp256k1.Fp(8) * u * u * u
                return EcJacobian(x3, y3, z3)
        else:
            u = u1 - u2
            v = v1 - v2
            w = u * u * z1 * z2 - v * v * v - secp256k1.Fp(2) * v * v * x1 * z2
            x3 = v * w
            y3 = u * (v * v * x1 * z2 - w) - v * v * v * y1 * z2
            z3 = v * v * v * z1 * z2
            return EcJacobian(x3, y3, z3)

    def __mul__(self, k):
        n = k.x
        result = secp256k1.I
        addend = self
        while n:
            b = n & 1
            if b == 1:
                result += addend
            addend = addend + addend
            n = n >> 1
        return result


p = secp256k1.G * secp256k1.Fr(42)
q = secp256k1.G * secp256k1.Fr(24)
assert p + q == secp256k1.G * secp256k1.Fr(66)
assert (EcJacobian.encode(p) + EcJacobian.encode(q)).decode() == secp256k1.G * secp256k1.Fr(66)

assert p + p == secp256k1.G * secp256k1.Fr(84)
assert (EcJacobian.encode(p) + EcJacobian.encode(p)).decode() == secp256k1.G * secp256k1.Fr(84)

r = secp256k1.Ec(p.x, -p.y)
assert p + r == secp256k1.I
assert (EcJacobian.encode(p) + EcJacobian.encode(r)).decode() == secp256k1.I

assert p + secp256k1.I == p
assert (EcJacobian.encode(p) + EcJacobian.encode(secp256k1.I)).decode() == p

assert p * secp256k1.Fr(42) == secp256k1.G * secp256k1.Fr(1764)
assert (EcJacobian.encode(p) * secp256k1.Fr(42)).decode() == p * secp256k1.Fr(42)