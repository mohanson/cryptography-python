import numpy

# p(x) = 5x² - 2x + 4
px = [4, -2, 5]

# q(x) = 2x² - 5x + 2
qx = [2, -5, 2]

r = numpy.polynomial.polynomial.polyadd(px, qx)
print('polyadd', r)  # 7x² - 7x + 6
r = numpy.polynomial.polynomial.polysub(px, qx)
print('polysub', r)  # 3x² + 3x + 2
r = numpy.polynomial.polynomial.polymul(px, qx)
print('polymul', r)  # 10x⁴ - 29x³ + 28x² - 24x + 8
quo, rem = numpy.polynomial.polynomial.polydiv(px, qx)
print('polyquo', quo)  # 2.5
print('polyrem', rem)  # 10.5x - 1
