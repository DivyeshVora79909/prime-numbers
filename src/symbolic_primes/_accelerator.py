"""Dense signed rows exchanged with the optional Rust extension.

One byte per ordinary base-two position, least significant first. Digits
+1 and -1 are encoded as 1 and 255. Towers and sparse gaps stay symbolic.
"""

try:
    from . import _native as native
except ImportError:
    native = None


MAX_DENSE_SPAN = 1 << 23


def ordinary_terms(value):
    """Canonical signed positions of a nonnegative, ordinary exponent value."""
    # value = (3*value - value)/2. Only differing bits survive subtraction.
    occupied = ((3 * value) ^ value) >> 1
    negative = value >> 1
    while occupied:
        bit = occupied & -occupied
        yield bit.bit_length() - 1, -1 if negative & bit else 1
        occupied ^= bit


def dense_row(value, *, min_terms=1, offset=0):
    terms = value.terms
    if len(terms) < min_terms or not terms:
        return None
    highest = terms[-1][0]._small
    if highest is None:
        return None
    span = highest - offset + 1
    if span <= 0 or span > min(MAX_DENSE_SPAN, max(256, 16 * len(terms))):
        return None
    row = bytearray(span)
    for exponent, sign in terms:
        position = exponent._small
        if position is None or position < offset:
            return None
        row[position - offset] = sign & 255
    return bytes(row)
