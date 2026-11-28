"""Bounded arithmetic on signed rows, independently of the exponent origin."""

from . import _accelerator

try:
    from numpy import convolve, int64, zeros
except ImportError:
    convolve = None


_NUMPY_SPAN_LIMIT = 32768
_NUMPY_WORK_LIMIT = 8_000_000


def combine_if_dense(circuit, a, b, sign):
    native = _accelerator.native
    if native is None or len(a.terms) + len(b.terms) < 64:
        return None
    start = min(a.terms[0][0]._small, b.terms[0][0]._small)
    for terms in (a.terms, b.terms):
        span = terms[-1][0]._small - start + 1
        if span > min(native.MAX_MULTIPLY_SPAN, max(256, 16 * len(terms))):
            return None
    left = _accelerator.dense_row(a, offset=start)
    right = _accelerator.dense_row(b, offset=start)
    result = native.add_rows(left, right, sign < 0)
    # Most output digits refer to one of the inputs' existing exponent nodes.
    # Reuse those nodes directly; create definitions only for new positions.
    positions = {e._small: e for marking in (a, b) for e, _ in marking.terms}
    edges = []
    for i, digit in enumerate(result, start):
        if digit:
            exponent = positions.get(i)
            if exponent is None:
                exponent = circuit.integer(i)
            edges.append((exponent, 1 if digit == 1 else -1))
    return circuit._intern(tuple(edges))


def dense_groups(groups):
    """Encode reusable local rows inside groups with symbolic high exponents."""
    rows = {}
    native = _accelerator.native
    if native is None:
        return rows
    for high, offsets in groups.items():
        start = offsets[0][0]
        span = offsets[-1][0] - start + 1
        if len(offsets) < 16 or span > min(native.MAX_MULTIPLY_SPAN, 6 * len(offsets)):
            continue
        row = bytearray(span)
        for offset, sign in offsets:
            row[offset - start] = sign & 255
        rows[high] = start, bytes(row)
    return rows


def accumulate_dense(left, right, coefficients):
    """Add a normalized signed product to one shared high-part accumulator."""
    left_base, left_row = left
    right_base, right_row = right
    native = _accelerator.native
    if len(left_row) + len(right_row) - 1 > native.MAX_MULTIPLY_SPAN:
        return False
    base = left_base + right_base
    for i, digit in enumerate(native.multiply_rows(left_row, right_row)):
        if digit:
            offset = base + i
            coefficients[offset] = coefficients.get(offset, 0) + (1 if digit == 1 else -1)
    return True


def multiply_if_dense(circuit, a, b):
    """Convolve bounded coefficient arrays, independent of the base exponent."""
    pairs = len(a.terms) * len(b.terms)
    native = _accelerator.native
    if pairs < 128 or (native is None and (convolve is None or pairs < 2048)):
        return None
    start_a, start_b = a.terms[0][0]._small, b.terms[0][0]._small
    span_a = a.terms[-1][0]._small - start_a + 1
    span_b = b.terms[-1][0]._small - start_b + 1
    if span_a > 6 * len(a.terms) or span_b > 6 * len(b.terms):
        return None
    base = start_a + start_b
    if native is not None and span_a + span_b - 1 <= native.MAX_MULTIPLY_SPAN:
        left = _accelerator.dense_row(a, offset=start_a)
        right = _accelerator.dense_row(b, offset=start_b)
        if left is None or right is None:
            return None
        product = native.multiply_rows(left, right)
        return circuit._from_small_positions(
            (base + i, 1 if digit == 1 else -1)
            for i, digit in enumerate(product) if digit
        )

    if (convolve is None or span_a + span_b > _NUMPY_SPAN_LIMIT
            or span_a * span_b > _NUMPY_WORK_LIMIT):
        return None
    left, right = zeros(span_a, dtype=int64), zeros(span_b, dtype=int64)
    for exponent, sign in a.terms:
        left[exponent._small - start_a] = sign
    for exponent, sign in b.terms:
        right[exponent._small - start_b] = sign
    # Each exact integer coefficient has magnitude at most the shorter span.
    coefficients = convolve(left, right).tolist()
    coefficients.append(0)
    terms = []
    i = 0
    while i < len(coefficients) - 1:
        coefficient = coefficients[i]
        if coefficient:
            digit = 0
            if coefficient % 2:
                digit = 2 - (coefficient + 2 * coefficients[i + 1]) % 4
                terms.append((base + i, digit))
            coefficients[i + 1] += (coefficient - digit) // 2
            if i + 2 == len(coefficients) and coefficients[-1]:
                coefficients.append(0)
        i += 1
    return circuit._from_small_positions(terms)
