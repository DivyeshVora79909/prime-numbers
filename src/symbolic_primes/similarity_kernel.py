"""Exact signed convolution using shared high parts of recursive exponents."""

from .symbolic import Circuit, Signed, _naf
from ._dense_multiply import accumulate_dense, dense_groups

_OFFSET_BITS = 64
_SAFE_OFFSET = 1 << (_OFFSET_BITS - 2)
_MAX_MASK_COORDINATES = 2048


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
    if a.terms[-1][0]._small is not None and b.terms[-1][0]._small is not None:
        return circuit.multiply(a, b)
    left = _exponent_groups(circuit, a)
    right = _exponent_groups(circuit, b)
    left_rows, right_rows = dense_groups(left), dense_groups(right)
    coordinates = _high_coordinates(circuit, left.keys() | right.keys())
    high_sums = {}
    blocks = {}
    for high_left, left_offsets in left.items():
        for high_right, right_offsets in right.items():
            high = _sum_high(circuit, high_left, high_right, coordinates, high_sums)
            coefficients = blocks.setdefault(high, {})
            if (high_left in left_rows and high_right in right_rows
                    and accumulate_dense(left_rows[high_left], right_rows[high_right], coefficients)):
                continue
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
