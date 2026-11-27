"""Signed hereditary base-2 arithmetic and compact exact quotients."""

from dataclasses import dataclass
from functools import cache
from heapq import heapify, heappop, heappush
from weakref import WeakValueDictionary


@dataclass(frozen=True, eq=False, slots=True, weakref_slot=True, repr=False)
class Signed:
    """An interned marking, ordered from smallest to largest exponent."""
    terms: tuple[tuple["Signed", int], ...]
    _small: int | None
    _owner: object

    def __lt__(self, other: "Signed") -> bool:
        return _compare(self, other) < 0


def _compare(a: Signed, b: Signed) -> int:
    if a is b:
        return 0
    if a._small is not None and b._small is not None:
        return (a._small > b._small) - (a._small < b._small)
    parents = []
    i, j = len(a.terms) - 1, len(b.terms) - 1
    while True:
        while i >= 0 and j >= 0:
            x, sx = a.terms[i]
            y, sy = b.terms[j]
            if x is not y:
                if x._small is None or y._small is None:
                    parents.append((a, b, i, j))
                    a, b = x, y
                    i, j = len(x.terms) - 1, len(y.terms) - 1
                    continue
                if x._small != y._small:
                    order = sx if x._small > y._small else -sy
                    break
            if sx != sy:
                order = sx - sy
                break
            i -= 1
            j -= 1
        else:
            order = 0
            if i >= 0:
                order = a.terms[i][1]
            elif j >= 0:
                order = -b.terms[j][1]

        # Propagate the exponent comparison back to its signed parent terms.
        # Equal terms resume the parent's comparison at the next lower digit.
        while parents:
            a, b, i, j = parents.pop()
            sx, sy = a.terms[i][1], b.terms[j][1]
            if order:
                order = sx if order > 0 else -sy
            elif sx != sy:
                order = sx - sy
            else:
                i -= 1
                j -= 1
                break
        else:
            return order


def _naf(coefficients, successor):
    """Consume sparse coefficients and emit nonadjacent signed positions."""
    queue = list(coefficients)
    heapify(queue)
    while queue:
        exponent = heappop(queue)
        coefficient = coefficients.pop(exponent)
        if not coefficient:
            continue
        if not queue and coefficient in (-1, 1):
            yield exponent, coefficient
            return
        following = successor(exponent)
        digit = 0
        if coefficient % 2:
            next_coefficient = coefficients.get(following, 0)
            digit = 2 - (coefficient + 2 * next_coefficient) % 4
            yield exponent, digit
        carry = (coefficient - digit) // 2
        if carry:
            if following not in coefficients:
                heappush(queue, following)
            coefficients[following] = coefficients.get(following, 0) + carry


def _varint(data: bytearray, n: int) -> None:
    while n >= 128:
        data.append((n & 127) | 128)
        n >>= 7
    data.append(n)


def _encode_rows(rows) -> bytes:
    data = bytearray(b"H2N1")
    _varint(data, len(rows))
    for row in rows:
        _varint(data, len(row))
        for child, digit in row:
            _varint(data, 2 * child + (digit < 0))
    return bytes(data)


def _encode_graph_rows(root, row_for):
    """Number each shared node after its children, in signed-row order."""
    indices = {}
    rows = []
    stack = [(root, False)]
    while stack:
        node, visited = stack.pop()
        if node in indices:
            continue
        terms = row_for(node)
        if visited:
            indices[node] = len(rows)
            rows.append([(indices[exponent], digit) for exponent, digit in terms])
        else:
            stack.append((node, True))
            stack.extend((exponent, False) for exponent, _ in reversed(terms))
    return rows


class Circuit:
    """An arena for canonical markings. Compare roots within the same arena."""

    def __init__(self):
        self._owner = object()
        self._nodes = WeakValueDictionary()
        self.zero = self._intern(())
        self.one = self._intern(((self.zero, 1),))

    def _check(self, *values: Signed) -> None:
        for value in values:
            if not isinstance(value, Signed) or value._owner is not self._owner:
                raise ValueError("expected a marking from this circuit")

    def _intern(self, terms) -> Signed:
        node = self._nodes.get(terms)
        if node is None:
            highest = terms[-1][0]._small if terms else 0
            small = None
            if highest is not None and highest <= 62:
                small = sum(digit * (1 << exponent._small) for exponent, digit in terms)
            node = Signed(terms, small, self._owner)
            self._nodes[terms] = node
        return node

    def integer(self, n: int) -> Signed:
        """Import an ordinary input integer, typically a small ID or position."""
        if not isinstance(n, int):
            raise TypeError("expected an int")
        terms = []
        exponent = 0
        while n:
            gap = (n & -n).bit_length() - 1
            n >>= gap
            exponent += gap
            digit = 2 - (n % 4)
            terms.append((self.integer(exponent), digit))
            n -= digit
        return self._intern(tuple(terms))

    def power(self, exponent: Signed) -> Signed:
        self._check(exponent)
        if exponent.terms and exponent.terms[-1][1] < 0:
            raise ValueError("negative power exponent")
        return self._intern(((exponent, 1),))

    def compare(self, a: Signed, b: Signed) -> int:
        self._check(a, b)
        return _compare(a, b)

    def _from_small_positions(self, terms) -> Signed:
        return self._intern(tuple((self.integer(exponent), digit) for exponent, digit in terms))

    def _from_coefficients(self, coefficients, small: bool) -> Signed:
        if small:
            return self._from_small_positions(_naf(coefficients, lambda e: e + 1))
        return self._intern(tuple(_naf(coefficients, lambda e: self.add(e, self.one))))

    def _combine(self, a: Signed, b: Signed, sign: int) -> Signed:
        small = all(e._small is not None for marking in (a, b) for e, _ in marking.terms)
        coefficients = {}
        for marking, multiplier in ((a, 1), (b, sign)):
            for exponent, digit in marking.terms:
                position = exponent._small if small else exponent
                coefficients[position] = coefficients.get(position, 0) + multiplier * digit
        return self._from_coefficients(coefficients, small)

    def add(self, a: Signed, b: Signed) -> Signed:
        self._check(a, b)
        if a._small is not None and b._small is not None:
            return self.integer(a._small + b._small)
        if a is self.one or b is self.one:
            return self._successor(b if a is self.one else a)
        return self._combine(a, b, 1)

    def _successor(self, value: Signed) -> Signed:
        """Carry a unit through the affected low digits, preserving the suffix."""
        terms = value.terms
        index = bit = 0
        carry = 1
        position = self.zero
        prefix = []
        while carry:
            coefficient = carry
            if index < len(terms) and terms[index][0] is position:
                coefficient += terms[index][1]
                index += 1
            bit += 1
            following = self.integer(bit)
            digit = 0
            if coefficient % 2:
                next_digit = 0
                if index < len(terms) and terms[index][0] is following:
                    next_digit = terms[index][1]
                digit = 2 - (coefficient + 2 * next_digit) % 4
                prefix.append((position, digit))
            carry = (coefficient - digit) // 2
            position = following
        return self._intern(tuple(prefix) + terms[index:])

    def subtract(self, a: Signed, b: Signed) -> Signed:
        self._check(a, b)
        if a is b:
            return self.zero
        if a._small is not None and b._small is not None:
            return self.integer(a._small - b._small)
        return self._combine(a, b, -1)

    def multiply(self, a: "Expression", b: "Expression") -> "Expression":
        if not isinstance(a, Signed) or not isinstance(b, Signed):
            left, left_divisor = _rational_parts(self, a)
            right, right_divisor = _rational_parts(self, b)
            if left_divisor == 1 and left is self.integer(right_divisor):
                return right
            if right_divisor == 1 and right is self.integer(left_divisor):
                return left
            return _pack_exact(self, self.multiply(left, right), left_divisor * right_divisor)
        self._check(a, b)
        if not a.terms or not b.terms:
            return self.zero
        small = all(e._small is not None for marking in (a, b) for e, _ in marking.terms)
        coefficients = {}
        if small:
            left = tuple((e._small, d) for e, d in a.terms)
            right = tuple((e._small, d) for e, d in b.terms)
            for x, sx in left:
                for y, sy in right:
                    exponent = x + y
                    coefficients[exponent] = coefficients.get(exponent, 0) + sx * sy
        else:
            for x, sx in a.terms:
                for y, sy in b.terms:
                    exponent = self.add(x, y)
                    coefficients[exponent] = coefficients.get(exponent, 0) + sx * sy
        return self._from_coefficients(coefficients, small)

    def divide_power(self, a: Signed, exponent: Signed) -> Signed:
        """Exactly divide a by 2**exponent, then preserve the NAF marking."""
        self._check(a, exponent)
        if self.compare(exponent, self.zero) < 0:
            raise ValueError("negative divisor exponent")
        if a.terms and self.compare(a.terms[0][0], exponent) < 0:
            raise ValueError("division is not exact")
        return self._intern(tuple((self.subtract(e, exponent), d) for e, d in a.terms))

    def bit_length(self, value: Signed) -> Signed:
        """Return the bit length of the absolute value as a recursive marking."""
        self._check(value)
        if not value.terms:
            return self.zero
        highest, sign = value.terms[-1]
        if len(value.terms) > 1 and value.terms[-2][1] != sign:
            return highest
        return self.add(highest, self.one)

    def _encoding(self, value: Signed):
        self._check(value)
        return _encode_graph_rows(value, lambda node: node.terms)

    def encode(self, value: Signed) -> bytes:
        """Encode a pure DAG as bottom-up varint rows. The last row is the root."""
        return _encode_rows(self._encoding(value))

    def stats(self, value: Signed) -> dict:
        rows = self._encoding(value)
        return {"terms": len(value.terms), "nodes": len(rows),
                "edges": sum(map(len, rows)), "bytes": len(_encode_rows(rows))}


@dataclass(frozen=True, slots=True)
class Geometric:
    """The quotient (2**(2*count) - 1)/3, with a recursive term count."""
    count: Signed


@dataclass(frozen=True, slots=True)
class ExactQuotient:
    """A checked quotient with an explicit positive denominator."""
    numerator: Signed
    denominator: int


Expression = Signed | Geometric | ExactQuotient


def _rational_parts(circuit: Circuit, value: Expression) -> tuple[Signed, int]:
    if isinstance(value, Geometric):
        circuit._check(value.count)
        doubled = circuit.add(value.count, value.count)
        return circuit.subtract(circuit.power(doubled), circuit.one), 3
    if isinstance(value, ExactQuotient):
        circuit._check(value.numerator)
        return value.numerator, value.denominator
    circuit._check(value)
    return value, 1


def _bounded_integer(value: Signed, bits: int) -> int:
    if value.terms:
        highest = value.terms[-1][0]._small
        if highest is None or highest > bits:
            raise ValueError(f"marking exceeds {bits}-bit integer import limit")
    result = sum(digit * (1 << exponent._small) for exponent, digit in value.terms)
    if result.bit_length() > bits:
        raise ValueError(f"marking exceeds {bits}-bit integer import limit")
    return result


def _pack_exact(circuit: Circuit, numerator: Signed, denominator: int) -> Expression:
    if denominator == 1 or not numerator.terms:
        return numerator
    if (denominator == 3 and len(numerator.terms) == 2
            and numerator.terms[0] == (circuit.zero, -1)
            and numerator.terms[1][1] == 1):
        exponent = numerator.terms[1][0]
        if exponent.terms[0][0] is not circuit.zero:
            return Geometric(circuit.divide_power(exponent, circuit.one))
    return ExactQuotient(numerator, denominator)


def exact_divide(circuit: Circuit, numerator: Expression, denominator: int | Signed,
                 *, periods=None) -> Expression:
    """Return a pure shift, a geometric sum, or a certified exact quotient."""
    value, existing = _rational_parts(circuit, numerator)
    exponent = None
    if isinstance(denominator, Signed):
        circuit._check(denominator)
        if len(denominator.terms) == 1 and denominator.terms[0][1] == 1:
            exponent = denominator.terms[0][0]
        else:
            denominator = _bounded_integer(denominator, 512)
    if exponent is None:
        if not isinstance(denominator, int) or denominator <= 0:
            raise ValueError("denominator must be a positive int or marking")
        if denominator & (denominator - 1) == 0:
            exponent = circuit.integer(denominator.bit_length() - 1)
    if exponent is not None:
        shifted = circuit.divide_power(value, exponent)
        if shifted is value:
            return numerator
        if existing != 1 and residue(circuit, shifted, existing, periods=periods):
            raise ValueError("division is not exact")
        return _pack_exact(circuit, shifted, existing)
    combined = existing * denominator
    result = _pack_exact(circuit, value, combined)
    if isinstance(result, ExactQuotient) and residue(circuit, value, combined, periods=periods):
        raise ValueError("division is not exact")
    return result


def _totient(n: int) -> int:
    remaining, phi = n, n
    p = 3
    while p * p <= remaining:
        if remaining % p == 0:
            phi -= phi // p
            while remaining % p == 0:
                remaining //= p
        p += 2
    if remaining > 1:
        phi -= phi // remaining
    return phi


def residue(circuit: Circuit, value, modulus: int, *, periods=None) -> int:
    """Evaluate modulo a positive modulus with no large integer expansion.

    For a non-small exponent, automatic period discovery is limited to
    odd moduli below 2**32. Supply a valid period of 2 for larger odd parts
    if needed, using a mapping from odd modulus to period.
    """
    if not isinstance(modulus, int) or modulus <= 0:
        raise ValueError("modulus must be a positive int")
    periods = {} if periods is None else periods

    @cache
    def period(odd):
        if odd in periods:
            bound = periods[odd]
            if bound <= 0 or pow(2, bound, odd) != 1:
                raise ValueError("invalid supplied period for 2")
            return bound
        if odd >= 1 << 32:
            raise ValueError("recursive exponent needs a supplied period for this odd modulus")
        return _totient(odd)

    @cache
    def power_mod(exponent, m):
        if exponent._small is not None:
            return pow(2, exponent._small, m)
        twos = (m & -m).bit_length() - 1
        two_power = 1 << twos
        odd = m >> twos
        mod_two = 0
        if circuit.compare(exponent, circuit.integer(twos)) < 0:
            mod_two = 1 << _bounded_integer(exponent, twos.bit_length() + 1)
        if odd == 1:
            return mod_two
        mod_odd = pow(2, evaluate(exponent, period(odd)), odd)
        if not twos:
            return mod_odd
        return mod_two + two_power * ((mod_odd - mod_two) * pow(two_power, -1, odd) % odd)

    @cache
    def evaluate(node, m):
        if m == 1:
            return 0
        if node._small is not None:
            return node._small % m
        return sum(digit * power_mod(exponent, m) for exponent, digit in node.terms) % m

    if isinstance(value, Geometric):
        circuit._check(value.count)
        exponent = circuit.add(value.count, value.count)
        return (power_mod(exponent, 3 * modulus) - 1) // 3
    if isinstance(value, ExactQuotient):
        circuit._check(value.numerator)
        return evaluate(value.numerator, modulus * value.denominator) // value.denominator
    circuit._check(value)
    return evaluate(value, modulus)
