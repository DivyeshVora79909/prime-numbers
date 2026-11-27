"""Exact signed convolution using shared high parts of recursive exponents."""

from .symbolic import Circuit, Signed, _naf

try:
    from numpy import convolve, int64, zeros
except ImportError:
    convolve = None


_OFFSET_BITS = 64
_SAFE_OFFSET = 1 << (_OFFSET_BITS - 2)
_MAX_MASK_COORDINATES = 2048
_MAX_DENSE_SPAN = 32768
_MAX_DENSE_WORK = 8_000_000


def _dense_small(circuit, a, b):
    """Convolve bounded coefficient arrays, independent of the base exponent."""
    if convolve is None or len(a.terms) * len(b.terms) < 2048:
        return None
    start_a, start_b = a.terms[0][0]._small, b.terms[0][0]._small
    span_a = a.terms[-1][0]._small - start_a + 1
    span_b = b.terms[-1][0]._small - start_b + 1
    if (span_a + span_b > _MAX_DENSE_SPAN or span_a * span_b > _MAX_DENSE_WORK
            or span_a > 6 * len(a.terms) or span_b > 6 * len(b.terms)):
        return None
    left, right = zeros(span_a, dtype=int64), zeros(span_b, dtype=int64)
    for exponent, sign in a.terms:
        left[exponent._small - start_a] = sign
    for exponent, sign in b.terms:
        right[exponent._small - start_b] = sign
    # Each exact integer coefficient has magnitude at most the shorter span.
    coefficients = convolve(left, right).tolist()
    coefficients.append(0)
    base = start_a + start_b
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


def _exponent_groups(circuit, value):
    groups = {}
    for exponent, digit in value.terms:
        offset = 0
        split = 0
        for position, sign in exponent.terms:
            bit = position._small
            if bit is None or bit >= _OFFSET_BITS:
                break
            offset += sign << bit
            split += 1
        high = circuit._intern(exponent.terms[split:])
        groups.setdefault(high, []).append((offset, digit))
    return groups


def _high_coordinates(circuit, groups):
    occupied = {position for high in groups for position, _ in high.terms}
    if len(occupied) > _MAX_MASK_COORDINATES:
        return None
    positions = sorted(occupied)
    ranks = {position: rank for rank, position in enumerate(positions)}
    adjacent = 0
    for rank, (first, second) in enumerate(zip(positions, positions[1:])):
        if circuit.add(first, circuit.one) is second:
            adjacent |= 1 << rank
    signatures = {}
    for high in groups:
        positive = negative = 0
        for position, digit in high.terms:
            marker = 1 << ranks[position]
            if digit > 0:
                positive |= marker
            else:
                negative |= marker
        signatures[high] = positive, negative
    return positions, adjacent, signatures


def _sum_high(circuit, left, right, coordinates, cache):
    if not left.terms:
        return right
    if not right.terms:
        return left
    if coordinates is None:
        return circuit.add(left, right)
    positions, adjacent, signatures = coordinates
    left_positive, left_negative = signatures[left]
    right_positive, right_negative = signatures[right]
    if left_positive & right_positive or left_negative & right_negative:
        return circuit.add(left, right)
    opposite = (left_positive & right_negative) | (left_negative & right_positive)
    positive = (left_positive | right_positive) ^ opposite
    negative = (left_negative | right_negative) ^ opposite
    occupied = positive | negative
    if occupied & (occupied >> 1) & adjacent:
        return circuit.add(left, right)
    signature = positive, negative
    result = cache.get(signature)
    if result is None:
        terms = []
        while occupied:
            bit = occupied & -occupied
            terms.append((positions[bit.bit_length() - 1], 1 if bit & positive else -1))
            occupied ^= bit
        result = circuit._intern(tuple(terms))
        cache[signature] = result
    return result


def multiply_fast(circuit: Circuit, a: Signed, b: Signed) -> Signed:
    """Multiply without expanding either number or its large exponents.

    Each exponent is split into a multiple of 2**64 and a signed local offset.
    High-part masks index occupied coordinates of the input DAG, independent
    of the numerical sizes of those coordinates.
    """
    circuit._check(a, b)
    if not a.terms or not b.terms:
        return circuit.zero
    if all(exponent._small is not None for value in (a, b) for exponent, _ in value.terms):
        dense = _dense_small(circuit, a, b)
        return circuit.multiply(a, b) if dense is None else dense
    left = _exponent_groups(circuit, a)
    right = _exponent_groups(circuit, b)
    coordinates = _high_coordinates(circuit, left.keys() | right.keys())
    high_sums = {}
    blocks = {}
    for high_left, left_offsets in left.items():
        for high_right, right_offsets in right.items():
            high = _sum_high(circuit, high_left, high_right, coordinates, high_sums)
            coefficients = blocks.setdefault(high, {})
            for x, sx in left_offsets:
                for y, sy in right_offsets:
                    offset = x + y
                    coefficients[offset] = coefficients.get(offset, 0) + sx * sy

    normalized = {}
    separated = True
    for high, coefficients in blocks.items():
        local = tuple(_naf(coefficients, lambda offset: offset + 1))
        if not local:
            continue
        normalized[high] = local
        if local[0][0] <= -_SAFE_OFFSET or local[-1][0] >= _SAFE_OFFSET:
            separated = False

    offsets = {}
    if separated:
        # The normalized local marking and its high marking have a blank
        # coordinate between them. Adjacent blocks cannot touch or borrow.
        output = []
        for high in sorted(normalized):
            for offset, digit in normalized[high]:
                local = offsets.get(offset)
                if local is None:
                    local = circuit.integer(offset)
                    offsets[offset] = local
                exponent = circuit._intern(local.terms + high.terms)
                output.append((exponent, digit))
        return circuit._intern(tuple(output))

    # Large offsets and carries can cross a block boundary. The general
    # normalizer sees the combined absolute positions in this case.
    coefficients = {}
    for high, local_terms in normalized.items():
        for offset, digit in local_terms:
            local = offsets.get(offset)
            if local is None:
                local = circuit.integer(offset)
                offsets[offset] = local
            exponent = circuit.add(high, local)
            coefficients[exponent] = coefficients.get(exponent, 0) + digit
    return circuit._from_coefficients(coefficients, False)
