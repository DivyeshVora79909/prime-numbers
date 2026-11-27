"""Output-sensitive exact division of signed hereditary base-2 markings."""

from heapq import merge
from operator import itemgetter

from .symbolic import Circuit, Signed


class QuotientWorkLimit(RuntimeError):
    """The requested quotient digit budget was exhausted."""


def _negate(circuit, value):
    return circuit._intern(tuple((exponent, -digit) for exponent, digit in value.terms))


def _shift(circuit, value, exponent):
    if exponent is circuit.zero:
        return value
    return circuit._intern(tuple(
        (circuit.add(position, exponent), digit) for position, digit in value.terms
    ))


def _subtract(circuit, left, right):
    if left.terms[-1][0]._small is not None and right.terms[-1][0]._small is not None:
        return circuit.subtract(left, right)
    positive = ((exponent, digit, False) for exponent, digit in left.terms)
    negative = ((exponent, -digit, True) for exponent, digit in right.terms)
    remaining = []
    for position, digit, source in merge(positive, negative, key=itemgetter(0)):
        if remaining:
            previous, previous_digit, previous_source = remaining[-1]
            if previous is position:
                if previous_digit == digit:
                    return circuit.subtract(left, right)
                remaining.pop()
                continue
            if source != previous_source and circuit.add(previous, circuit.one) is position:
                return circuit.subtract(left, right)
        remaining.append((position, digit, source))
    return circuit._intern(tuple((position, digit) for position, digit, _ in remaining))


def _above_threshold(circuit, remainder, shifted, offset):
    if remainder._small is not None and shifted._small is not None:
        return 3 * remainder._small > shifted._small << offset._small

    def descending(marking, shift, sign):
        for exponent, digit in reversed(marking.terms):
            position = exponent if shift is circuit.zero else circuit.add(exponent, shift)
            yield position, sign * digit

    streams = (
        descending(remainder, circuit.zero, 1),
        descending(remainder, circuit.one, 1),
        descending(shifted, offset, -1),
    )
    heads = [next(stream, None) for stream in streams]
    previous = None
    coefficient = 0
    while any(head is not None for head in heads):
        position = max(head[0] for head in heads if head is not None)
        if coefficient:
            # The three NAF tails below the processed prefix have
            # combined absolute value strictly below 4*2**position.
            if abs(coefficient) >= 2 or circuit.add(position, circuit.one) is not previous:
                return coefficient > 0
            coefficient *= 2
        for i, head in enumerate(heads):
            if head is not None and head[0] is position:
                coefficient += head[1]
                heads[i] = next(streams[i], None)
        previous = position
    return coefficient > 0


def divide_fast(
    circuit: Circuit, numerator: Signed, denominator: Signed, *, max_terms: int = 256
) -> Signed:
    """Return the pure canonical exact quotient, subject to a digit budget.

    The numerator and denominator are symbolic integers of either sign.
    Inexact division raises ValueError if encountered within the budget.
    QuotientWorkLimit leaves integrality undecided when the work runs out.
    """
    circuit._check(numerator, denominator)
    if not isinstance(max_terms, int) or max_terms < 0:
        raise ValueError("max_terms must be a nonnegative integer")
    if not denominator.terms:
        raise ValueError("division by zero")
    if not numerator.terms:
        return circuit.zero

    denominator_sign = denominator.terms[-1][1]
    if len(denominator.terms) == 1:
        quotient = circuit.divide_power(numerator, denominator.terms[0][0])
        if len(quotient.terms) > max_terms:
            raise QuotientWorkLimit(f"quotient exceeded {max_terms} signed digits")
        return quotient if denominator_sign > 0 else _negate(circuit, quotient)

    numerator_sign = numerator.terms[-1][1]
    remainder = numerator if numerator_sign > 0 else _negate(circuit, numerator)
    divisor = denominator if denominator_sign > 0 else _negate(circuit, denominator)
    digit_sign = numerator_sign * denominator_sign
    divisor_length = circuit.bit_length(divisor)
    two = circuit.integer(2)
    digits = []

    while remainder is not circuit.zero:
        if circuit.compare(remainder, divisor) < 0:
            raise ValueError("division is not exact")
        if len(digits) == max_terms:
            raise QuotientWorkLimit(f"quotient exceeded {max_terms} signed digits")

        guess = circuit.subtract(circuit.bit_length(remainder), divisor_length)
        shifted = _shift(circuit, divisor, guess)
        order = circuit.compare(remainder, shifted)
        if order == 0:
            digits.append((guess, digit_sign))
            break

        # k is floor(log2(remainder / divisor)). The largest NAF ending
        # at exponent k is floor(2**(k+2) / 3).
        k = circuit.subtract(guess, circuit.one) if order < 0 else guess
        if _above_threshold(circuit, remainder, shifted, circuit.one if order < 0 else two):
            position = circuit.add(k, circuit.one)
            contribution = shifted if order < 0 else _shift(circuit, shifted, circuit.one)
            digits.append((position, digit_sign))
            remainder = _subtract(circuit, contribution, remainder)
            digit_sign = -digit_sign
        else:
            contribution = _shift(circuit, divisor, k) if order < 0 else shifted
            digits.append((k, digit_sign))
            remainder = _subtract(circuit, remainder, contribution)

    return circuit._intern(tuple(reversed(digits)))
