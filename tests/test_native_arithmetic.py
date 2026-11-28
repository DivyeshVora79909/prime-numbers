"""Regression cases for the signed-coordinate remainder and exact division kernels."""

from random import Random
import unittest
from unittest.mock import patch

from symbolic_primes import Adaptive, Circuit, divide_fast, residue
from symbolic_primes import _accelerator
from symbolic_primes.quotient_kernel import QuotientWorkLimit, _divide_dense


def outcome(operation):
    try:
        return operation()
    except (ValueError, QuotientWorkLimit) as exc:
        return type(exc)


@unittest.skipUnless(_accelerator.native is not None, "Rust extension not installed")
class NativeArithmeticTests(unittest.TestCase):
    def test_noncanonical_coefficients_and_negative_borrow_across_word_boundaries(self):
        rng = Random(23)
        c = Circuit()
        for span in (63, 64, 65, 128, 129, 257):
            rows = [
                b'\xff' + bytes(span - 2) + b'\1',
                b'\1' + bytes(span - 2) + b'\xff',
                bytes(rng.choice((0, 1, 255)) for _ in range(span)),
            ]
            for row in rows:
                value = sum((-1 if d == 255 else d) << i for i, d in enumerate(row))
                for m in (1, 9, 1 << 64, (1 << 131) - 13):
                    divisor = m.to_bytes((m.bit_length() + 7) // 8, "little")
                    got = _accelerator.native.dense_residue(row, divisor)
                    self.assertEqual(int.from_bytes(got, "little"), value % m)
                    quotient = _accelerator.native.dense_divide(row, divisor)
                    if value % m:
                        self.assertIsNone(quotient)
                    else:
                        expected = c.integer(value // m)
                        actual_terms = [(i, -1 if d == 255 else d)
                                        for i, d in enumerate(quotient) if d]
                        self.assertEqual(actual_terms, [(e._small, d) for e, d in expected.terms])

    def test_pure_quotient_budget_and_negative_signs(self):
        c = Circuit()
        q = Random(91).getrandbits(450)
        p = 4_909_305_867_477_490_416_066_453_498_166_360_224_791
        for numerator_sign, denominator_sign in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
            n, d = c.integer(numerator_sign * q * p), c.integer(denominator_sign * p)
            expected = c.integer(numerator_sign * denominator_sign * q)
            count = len(expected.terms)
            self.assertIs(divide_fast(c, n, d, max_terms=count), expected)
            with self.assertRaises(QuotientWorkLimit):
                divide_fast(c, n, d, max_terms=count-1)
            a = Adaptive(c)
            self.assertEqual(a.residue(a.view(n), p+6), numerator_sign * q * p % (p+6))
            self.assertEqual(residue(c, n, p), 0)

    def test_budget_boundary_matches_descending_division(self):
        c = Circuit()
        # Test every possible remainder, including exact and fractions on
        # opposite sides of the integer algorithm's signed-digit threshold.
        for divisor in (3, 5, 6, 7):
            d = c.integer(divisor)
            for value in range(1, 12 * divisor):
                n = c.integer(value)
                row = _accelerator.dense_row(n)
                for budget in range(4):
                    with patch.object(_accelerator, 'native', None):
                        expected = outcome(lambda: divide_fast(c, n, d, max_terms=budget))
                    actual = outcome(lambda: _divide_dense(row, divisor, budget))
                    if actual is None:
                        # Let the exact descending algorithm resolve a boundary.
                        actual = outcome(lambda: divide_fast(c, n, d, max_terms=budget))
                    elif isinstance(actual, bytes):
                        actual = c._from_small_positions(
                            (i, -1 if bit == 255 else 1) for i, bit in enumerate(actual) if bit
                        )
                    self.assertIs(actual, expected, (value, divisor, budget))
