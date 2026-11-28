"""Whole-file candidate and certificate for large, dense signed roots."""

from . import _accelerator
from ._accelerator import MAX_DENSE_SPAN, ordinary_terms
from ._storage_bound import _naf_weights, _root_digits
from .symbolic import _varint


def _encode(row):
    """H2N1 postorder over ordinary exponent definitions and one signed root."""
    ids = {}
    body = bytearray()
    top = bytearray()

    def intern(e):
        if e not in ids:
            refs = [2 * intern(i) + (sign < 0) for i, sign in ordinary_terms(e)]
            ids[e] = len(ids)
            _varint(body, len(refs))
            for ref in refs:
                _varint(body, ref)
        return ids[e]

    _varint(top, len(row) - row.count(0))
    for e, digit in enumerate(row):
        if digit:
            _varint(top, 2 * intern(e) + (digit == 255))
    header = bytearray(b"H2N1")
    _varint(header, len(ids) + 1)
    return bytes(header + body + top)


def wide_candidate(row):
    if not row or len(row) > MAX_DENSE_SPAN:
        return None
    minimum = len(row) - row.count(0)
    if minimum < 8192 or len(row) > 16 * minimum:
        return None
    if _accelerator.native is not None:
        return _accelerator.native.pack_dense(row)

    # The reference backend uses the same DP and tie rule as Rust.
    sign = -1 if row[-1] == 255 else 1
    bits = bytearray(len(row) + 1)
    carry = 0
    for e, digit in enumerate(row):
        carry += sign * ((digit == 1) - (digit == 255))
        bits[e] = carry & 1
        carry >>= 1
    costs = bytes(4 + weight for weight in _naf_weights(len(row)))
    chosen = _root_digits(bits, costs)
    alternative = bytearray(len(bits))
    cost = 0
    for e, digit in chosen:
        alternative[e] = digit & 255
        cost += costs[e]

    # Above the trailing zeros, rewrite the odd prefix P of width w using
    # +/-1 in every position. Its positive mask B solves 2B - (2**w - 1) = P.
    k = len(bits).bit_length()
    start = next((i for i in range(k + 1) if alternative[i]), None)
    if start is not None:
        low = alternative[start:k + 1]
        prefix = sum(((d == 1) - (d == 255)) << i for i, d in enumerate(low))
        positive = (prefix + (1 << len(low)) - 1) // 2
        for i in range(len(low)):
            alternative[start + i] = 1 if positive >> i & 1 else 255
    if sign < 0:
        alternative = alternative.translate(bytes.maketrans(b"\x01\xff", b"\xff\x01"))
    # Each selected coordinate needs a root reference, an arity, and at
    # least its NAF weight in defining edges. Count unavoidable wide IDs.
    header = 4 + (minimum.bit_length() + 6) // 7 + ((minimum + 1).bit_length() + 6) // 7
    overhead = header - 64 - 8192 + max(0, minimum - (1 << 20))
    best = _encode(alternative)
    canonical_cost = sum(weight for digit, weight in zip(row, costs) if digit)
    if len(best) >= overhead + canonical_cost:
        best = min((_encode(row), best), key=len)
    return best, overhead + cost
