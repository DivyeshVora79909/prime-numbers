"""Regression checks for compact rewrites, arithmetic and serialization."""

import unittest
from random import Random

from symbolic_primes.adaptive import Adaptive
from symbolic_primes.quotient_kernel import QuotientWorkLimit, divide_fast
from symbolic_primes.similarity_kernel import multiply_fast
from symbolic_primes.symbolic import Circuit


class AdaptiveTests(unittest.TestCase):
    def setUp(self):
        self.c = Circuit()
        self.a = Adaptive(self.c)

    def test_signed_addition_and_subtraction_at_translated_positions(self):
        c = self.c
        rng = Random(725)
        x = rng.getrandbits(720)
        y = -(rng.getrandbits(413) << 8)

        def translate(value, offset):
            return c._intern(tuple((c.integer(e._small + offset), d)
                                   for e, d in c.integer(value).terms))

        for offset in (0, 1 << 48):
            a, b = translate(x, offset), translate(y, offset)
            self.assertIs(c.add(a, b), translate(x + y, offset))
            self.assertIs(c.subtract(b, a), translate(y - x, offset))
            self.assertIs(c.add(a, translate(-x, offset)), c.zero)
            self.assertIs(c.subtract(c.zero, a), translate(-x, offset))
            self.assertIs(c.subtract(c.add(a, b), a), b)

        # A distant ordinary digit and an unbounded tower take the sparse path.
        for exponent in (c.integer(1 << 48), c.power(c.power(c.integer(80)))):
            far = c.power(exponent)
            dense = c.integer(x)
            self.assertIs(c.subtract(c.add(dense, far), far), dense)

    def test_dense_and_recursive_multiplication(self):
        c = self.c
        rng = Random(1937)
        for _ in range(16):
            x = rng.randrange(-(1 << 350), 1 << 350)
            y = rng.randrange(-(1 << 210), 1 << 210)
            a, b = c.integer(x), c.integer(y)
            expected = c.integer(x * y)
            self.assertIs(c.multiply(a, b), expected)
            self.assertIs(multiply_fast(c, a, b), expected)
        exponent = c.power(c.power(c.integer(80)))
        second = c.power(c.add(exponent, c.integer(15)))
        left = c.zero
        right = c.zero
        for offset in (5, 16, 43, 700):
            left = c.add(left, c.power(c.add(exponent, c.integer(offset))))
            right = c.subtract(right, c.power(c.add(second, c.integer(offset))))
        self.assertIs(multiply_fast(c, left, right), c.multiply(left, right))

    def test_dense_offset_products_at_tower_exponents(self):
        c = self.c
        high = c.power(c.power(c.integer(80)))
        twice_high = c.add(high, high)

        def translated(base, sign):
            return c._intern(tuple(
                (c.add(base, c.integer(offset)), sign if offset % 5 else -sign)
                for offset in range(-40, 58, 3)
            ))

        # Two cross-products land in the same high group; their carries interact.
        left = c.add(translated(high, 1), translated(twice_high, -1))
        right = c.add(translated(high, -1), translated(twice_high, 1))
        self.assertIs(multiply_fast(c, left, right), c.multiply(left, right))

    def test_repeated_block_in_an_exponent(self):
        c, a = self.c, self.a
        block = c.integer(0b10100001001)
        stride = c.power(c.power(c.integer(19)))
        position = c.zero
        translations = c.zero
        for index in range(12):
            term = c.power(position)
            translations = c.subtract(translations, term) if index % 2 else c.add(translations, term)
            position = c.add(position, stride)
        marking = c.power(c.subtract(c.zero, multiply_fast(c, block, translations)))
        packed = a.pack(marking)
        self.assertLess(a.stats(packed)["bytes"], a.stats(a.view(marking))["bytes"])
        self.assertIs(a.canonicalize(packed), marking)
        new = Adaptive()
        serialized = a.encode(packed)
        restored = new.decode(serialized)
        self.assertEqual(new.encode(restored), serialized)
        for modulus in (16, 97, 192):
            self.assertEqual(new.residue(restored, modulus), a.residue(packed, modulus))

    def test_whole_dag_cost_compared_with_literal(self):
        c, a = self.c, self.a
        fixture = c.integer(11281552818663026899)
        result = a.pack(fixture)
        self.assertLessEqual(a.stats(result)["bytes"], a.stats(a.view(fixture))["bytes"])
        self.assertIs(a.canonicalize(result), fixture)

    def test_persistent_updates_and_cancellation(self):
        c, a = self.c, self.a
        factors = tuple(a.view(c.integer(n)) for n in range(3, 515))
        product = a.one
        oracle = 1
        modulus = 65537
        for i, factor in enumerate(factors):
            product = a.multiply(product, factor)
            oracle = oracle * (i + 3) % modulus
        self.assertLessEqual(product._height, 12)
        self.assertEqual(a.residue(product, modulus), oracle)
        without = a.divide(product, factors[183])
        self.assertEqual(a.residue(without, modulus), oracle * pow(186, -1, modulus) % modulus)
        self.assertEqual(a.residue(a.multiply(without, factors[183]), modulus), oracle)
        self.assertEqual(a.residue(product, modulus), oracle)
        serialized = a.encode(product)
        new = Adaptive()
        self.assertEqual(new.residue(new.decode(serialized), modulus), oracle)
        with self.assertRaises(ValueError):
            new.decode(serialized[:-1])

    def test_repeat_geometric_and_shift_rules(self):
        c, a = self.c, self.a
        factor = a.view(c.integer(1337))
        cubed = a.multiply(a.multiply(factor, factor), factor)
        self.assertEqual((cubed.kind, cubed.payload), ("repeat", 3))
        self.assertEqual(a.divide(cubed, factor).payload, 2)
        stride = c.power(c.integer(64))
        series_a = c.zero
        series_b = c.zero
        position = c.zero
        for index in range(12):
            term = c.power(position)
            series_a = c.add(series_a, term)
            if index < 4:
                series_b = c.add(series_b, term)
            position = c.add(position, stride)
        full, part = map(a.pack, (series_a, series_b))
        self.assertEqual((full.kind, part.kind), ("geo", "geo"))
        after_division = a.divide(full, part)
        self.assertIs(a.multiply(after_division, part), full)
        body = c.integer(29)
        larger = c.power(c.add(stride, c.integer(12)))
        smaller = c.power(stride)
        numerator = a.pack(multiply_fast(c, larger, body))
        denominator = a.pack(multiply_fast(c, smaller, body))
        self.assertIs(a.canonicalize(a.divide(numerator, denominator)), c.integer(4096))

    def test_division_without_product_nodes(self):
        c, a = self.c, self.a
        large_exponent = c.power(c.power(c.integer(72)))
        n1 = c.subtract(c.power(large_exponent), c.integer(29))
        n2 = c.add(c.power(c.add(large_exponent, c.integer(117))), c.integer(19))
        numerator = multiply_fast(c, n1, n2)
        self.assertIs(divide_fast(c, numerator, n1), n2)
        self.assertIs(a.canonicalize(a.divide(a.view(numerator), a.view(n1))), n2)
        with self.assertRaises(ValueError):
            divide_fast(c, c.add(numerator, c.one), n1)
        with self.assertRaises(QuotientWorkLimit):
            divide_fast(c, numerator, n1, max_terms=1)
        if (quotient := a.divide(a.view(c.integer(45)), 3, max_terms=0)).kind != "quotient":
            self.fail("expected a checked quotient node")
        self.assertEqual(a.residue(quotient, 3), 0)
        self.assertIs(a.multiply(quotient, a.view(c.integer(3))), a.view(c.integer(45)))


if __name__ == "__main__":
    unittest.main()
