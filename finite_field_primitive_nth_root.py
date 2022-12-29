import finite_field
import random


class Fp(finite_field.Fp):

    @classmethod
    def primitive_nth_root(cls, n):
        # https://crypto.stackexchange.com/questions/63614/finding-the-n-th-root-of-unity-in-a-finite-field
        assert (cls.p - 1) % n == 0
        while 1:
            x = cls(random.randint(1, cls.p - 1))
            g = x ** ((cls.p - 1) // n)
            if g ** (n // 2) != Fp(1):
                return g


if __name__ == '__main__':
    Fp.p = 1 + 407 * (1 << 119)
    assert Fp.primitive_nth_root(1 << 119) ** (1 << 119) == Fp(1)
    assert Fp.primitive_nth_root(1 << 111) ** (1 << 111) == Fp(1)
    Fp.p = 0
