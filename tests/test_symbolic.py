"""Bounded reference checks and symbolic tower regression checks."""

import unittest
from random import Random

from symbolic_primes.symbolic import Circuit, ExactQuotient, Geometric, exact_divide, residue


class CircuitTests(unittest.TestCase):
    def setUp(self):
        self.c = Circuit()

    def test_canonical_recoding(self):
        c = self.c
        value = c.zero
        for index in (8, 6, 5, 4, 3, 2, 0):
            value = c.add(value, c.power(c.integer(index)))
        self.assertIs(value, c.integer(381))
        self.assertEqual(
            value.terms,
            tuple((c.integer(index), sign) for index, sign in [(0, 1), (2, -1), (7, -1), (9, 1)]),
        )

    def test_bounded_random_reference(self):
        c = self.c
        rng = Random(100)
        for _ in range(500):
            left = rng.randrange(-(1 << 100), 1 << 100)
            right = rng.randrange(-(1 << 100), 1 << 100)
            a, b = c.integer(left), c.integer(right)
            self.assertIs(c.bit_length(a), c.integer(left.bit_length()))
            self.assertIs(c.add(a, b), c.integer(left + right))
            self.assertIs(c.subtract(a, b), c.integer(left - right))
            self.assertIs(c.multiply(a, b), c.integer(left * right))
            order = c.compare(a, b)
            self.assertEqual((order > 0) - (order < 0), (left > right) - (left < right))
            divisor = rng.randrange(1, 1 << 30)
            numerator = c.integer(left * divisor)
            quotient = exact_divide(c, numerator, divisor)
            self.assertEqual(residue(c, quotient, divisor), left % divisor)

    def test_mixed_depth_carry_and_borrow(self):
        c = self.c
        lower = c.power(c.integer(64))
        exponent = c.power(lower)
        twice_exponent = c.power(c.add(lower, c.one))
        x = c.power(exponent)
        self.assertIs(c.add(exponent, exponent), twice_exponent)
        self.assertIs(c.add(x, x), c.power(c.add(exponent, c.one)))
        self.assertIs(c.subtract(c.power(c.add(exponent, c.one)), x), x)
        self.assertIs(c.subtract(c.power(c.add(exponent, c.one)), c.one),
                      c.add(c.subtract(x, c.one), x))
        product = c.multiply(c.subtract(x, c.one), c.add(x, c.one))
        self.assertIs(product, c.subtract(c.power(twice_exponent), c.one))
        self.assertIs(c.bit_length(product), twice_exponent)
        self.assertIs(exact_divide(c, c.multiply(product, x), x), product)
        self.assertEqual(c.stats(product), {"terms": 2, "nodes": 9, "edges": 12, "bytes": 26})

    def test_tower_division_and_noninvertible_residues(self):
        c = self.c
        count = c.power(c.power(c.integer(64)))
        numerator = c.subtract(c.power(c.add(count, count)), c.one)
        quotient = exact_divide(c, numerator, 3)
        self.assertIsInstance(quotient, Geometric)
        self.assertIs(quotient.count, count)
        self.assertEqual(residue(c, quotient, 3), 1)
        self.assertEqual(residue(c, quotient, 16), 5)
        self.assertIs(c.multiply(quotient, c.integer(3)), numerator)
        for modulus in (97, 96, 576, 768):
            self.assertEqual(residue(c, quotient, modulus) * 3 % modulus,
                             residue(c, numerator, modulus))

    def test_cascading_division(self):
        c = self.c
        numerator = c.integer(-45)
        first = exact_divide(c, numerator, c.integer(3))
        self.assertIsInstance(first, ExactQuotient)
        second = exact_divide(c, first, 5)
        self.assertEqual(residue(c, second, 9), 6)
        self.assertIs(c.multiply(second, c.integer(15)), numerator)
        shifted = exact_divide(c, exact_divide(c, c.integer(24), 3), 2)
        self.assertEqual(residue(c, shifted, 13), 4)
        self.assertRaises(ValueError, exact_divide, c, first, 2)
        self.assertRaises(ValueError, exact_divide, c, numerator, 7)
        self.assertRaises(ValueError, exact_divide, c, numerator, 0)

    def test_large_modulus_with_supplied_period(self):
        c = self.c
        exponent = c.power(c.power(c.integer(64)))
        value = c.power(exponent)
        modulus = 4_294_967_311
        expected = pow(2, pow(2, 1 << 64, modulus - 1), modulus)
        self.assertRaises(ValueError, residue, c, value, modulus)
        self.assertEqual(residue(c, value, modulus, periods={modulus: modulus - 1}),
                         expected)
        self.assertRaises(ValueError, residue, c, value, modulus, periods={modulus: 1})


if __name__ == "__main__":
    unittest.main()
