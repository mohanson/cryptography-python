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
