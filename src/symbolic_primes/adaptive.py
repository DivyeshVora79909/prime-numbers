"""Costed recursive factor discovery and exact lazy arithmetic."""

from collections import defaultdict
from dataclasses import dataclass
from functools import cache
from itertools import pairwise, product
from weakref import WeakValueDictionary

from .adaptive_io import decode_records, encode, residue, stats
from .quotient_kernel import QuotientWorkLimit, divide_fast
from .similarity_kernel import multiply_fast
from .symbolic import Circuit, Signed, _bounded_integer


@dataclass(frozen=True, eq=False, slots=True, weakref_slot=True, repr=False)
class Node:
    kind: str
    children: tuple["Node", ...]
    payload: int | None
    _owner: object
    _small: int | None
    _wide: bool
    _height: int
    _marker: int
    _filter: int


class _BudgetExceeded(Exception):
    pass


def _known_small(kind, children, payload):
    values = tuple(child._small for child in children)
    if any(value is None for value in values):
        return None
    if kind == "naf" and all(0 <= exponent <= 62 for exponent in values):
        number = sum((-1 if payload >> i & 1 else 1) << exponent
                     for i, exponent in enumerate(values))
    elif kind == "sum":
        number = sum(values)
    elif kind == "neg":
        number = -values[0]
    elif kind == "shift" and 0 <= values[1] <= 62:
        number = values[0] << values[1]
    elif kind == "mul":
        number = values[0] * values[1]
    elif kind == "repeat" and (abs(values[0]) <= 1 or payload <= 62):
        number = pow(values[0], payload)
    elif kind == "geo" and values[0] == 0:
        number = payload
    elif kind == "geo" and 0 < values[0] and values[0] * (payload - 1) <= 62:
        number = ((1 << (values[0] * payload)) - 1) // ((1 << values[0]) - 1) if payload > 1 else 1
    elif kind == "quotient" and values[0] % payload == 0:
        number = values[0] // payload
    else:
        return None
    return number if number.bit_length() <= 63 else None


class Adaptive:
    """One arena for compressed markings, exact factors and geometric runs."""

    def __init__(self, circuit: Circuit | None = None):
        self.circuit = Circuit() if circuit is None else circuit
        self._owner = object()
        self._nodes = WeakValueDictionary()
        self._serial = 0
        self.zero = self._node("naf", payload=0, small=0)
        self.one = self._node("naf", (self.zero,), 0, small=1)

    def _node(self, kind, children=(), payload=None, *, small=None, wide=False):
        key = kind, children, payload
        node = self._nodes.get(key)
        if node is None:
            self._serial += 1
            serial = self._serial
            bit = (serial * 0x9E3779B1 ^ (serial >> 5)) & 511
            marker = 1 << bit
            mask = marker
            if kind in ("mul", "repeat"):
                for child in children:
                    mask |= child._filter
            height = 1 + max(child._height for child in children) if kind == "mul" else 0
            node = Node(kind, children, payload, self._owner, small, wide, height, marker, mask)
            self._nodes[key] = node
        else:
            # These facts describe the value and are independent of node identity.
            if small is not None:
                object.__setattr__(node, "_small", small)
            if wide:
                object.__setattr__(node, "_wide", True)
        return node

    def _check(self, *values):
        for value in values:
            if not isinstance(value, Node) or value._owner is not self._owner:
                raise ValueError("expected an expression from this adaptive arena")

    def _naf(self, exponents, signs, raw):
        mask = sum((sign < 0) << i for i, sign in enumerate(signs))
        return self._node("naf", tuple(exponents), mask,
                          small=raw._small, wide=raw._small is None)

    def view(self, marking: Signed) -> Node:
        """Use the extended encoding with pure signed markings only."""
        self.circuit._check(marking)

        @cache
        def visit(raw):
            return self._naf((visit(e) for e, _ in raw.terms),
                             (d for _, d in raw.terms), raw)

        return visit(marking)

    def _repeat(self, value, count):
        if count == 1:
            return value
        if value.kind == "repeat":
            return self._repeat(value.children[0], count * value.payload)
        return self._node("repeat", (value,), count)

    def _shift(self, body, exponent):
        if body is self.zero or exponent is self.zero:
            return body
        if body is self.one:
            return self._node("naf", (exponent,), 0)
        return self._node("shift", (body, exponent))

    def _neg(self, value):
        if value is self.zero:
            return value
        if value.kind == "neg":
            return value.children[0]
        small = None if value._small is None else -value._small
        return self._node("neg", (value,), small=small, wide=value._wide)

    def _pair(self, a, b):
        return self._node("mul", (a, b))

    def _balance(self, a, b):
        if a._height > b._height + 1:
            left, right = a.children
            if left._height >= right._height:
                return self._pair(left, self._pair(right, b))
            middle, far = right.children
            return self._pair(self._pair(left, middle), self._pair(far, b))
        if b._height > a._height + 1:
            left, right = b.children
            if right._height >= left._height:
                return self._pair(self._pair(a, left), right)
            far, middle = left.children
            return self._pair(self._pair(a, far), self._pair(middle, right))
        return self._pair(a, b)

    def _join(self, a, b):
        if a is self.one:
            return b
        if b is self.one:
            return a
        if a._height > b._height + 1:
            left, right = a.children
            return self._balance(left, self._join(right, b))
        if b._height > a._height + 1:
            left, right = b.children
            return self._balance(self._join(a, left), right)
        return self._pair(a, b)

    def pack(self, marking: Signed, *, rounds: int = 1) -> Node:
        """Choose the smallest encoding found by bottom-up exact rewrites."""
        c = self.circuit
        c._check(marking)
        if rounds < 0:
            raise ValueError("rewrite rounds must be nonnegative")
        zero = c.zero

        @cache
        def byte_cost(expr):
            return len(encode(expr))

        @cache
        def difference(a, b):
            return c.subtract(a, b)

        @cache
        def literal(raw):
            return self._naf((literal(e) for e, _ in raw.terms),
                             (d for _, d in raw.terms), raw)

        @cache
        def visit(raw, budget):
            terms = raw.terms
            plain = literal(raw)
            choices = tuple(visit(e, rounds) for e, _ in terms)
            best = plain
            n = len(terms)

            def consider(candidate):
                nonlocal best
                if byte_cost(candidate) < byte_cost(best):
                    best = candidate

            if choices != plain.children:
                digits = tuple(d for _, d in terms)
                consider(self._naf(choices, digits, raw))
                changed = [i for i in range(n) if choices[i] is not plain.children[i]]
                # Local choices can lose sharing in their parent. Price the
                # entire DAG for every combination in this bounded portfolio.
                if 1 < len(changed) <= 6:
                    for switches in product((False, True), repeat=len(changed)):
                        if any(switches) and not all(switches):
                            children = list(plain.children)
                            for i, switch in zip(changed, switches):
                                if switch:
                                    children[i] = choices[i]
                            consider(self._naf(children, digits, raw))

            if n < 2 or budget == 0:
                return best

            exponents, digits = zip(*terms)
            gaps = tuple(difference(b, a) for a, b in pairwise(exponents))
            anchor = exponents[0]

            if anchor is not zero:
                span = difference(exponents[-1], anchor)
                if anchor._small is not None or c.compare(anchor, span) > 0:
                    relative = c._intern(tuple((difference(e, anchor), d) for e, d in terms))
                    consider(self._shift(visit(relative, budget - 1), visit(anchor, rounds)))

            if n >= 3 and len(set(digits)) == 1 and len(set(gaps)) == 1:
                series = self._node("geo", (visit(gaps[0], rounds),), n)
                if digits[0] < 0:
                    series = self._neg(series)
                consider(self._shift(series, visit(anchor, rounds)))

            # Equal internal gap and relative-sign sequences certify translation.
            @cache
            def signature(i, j):
                orientation = digits[i]
                return gaps[i:j - 1], tuple(orientation * d for d in digits[i:j])

            examined = set()

            def build(groups):
                grouping = tuple(tuple(group) for group in groups if len(group) >= 2)
                if not grouping or grouping in examined:
                    return
                examined.add(grouping)
                covered = set()
                pieces = []
                for group in grouping:
                    first, end = group[0]
                    origin = exponents[first]
                    orientation = digits[first]
                    motif = c._intern(tuple((difference(e, origin), d * orientation)
                                            for e, d in terms[first:end]))
                    anchors = c._intern(tuple(terms[i] for i, _ in group))
                    pieces.append(self.multiply(visit(anchors, budget - 1), visit(motif, budget - 1)))
                    for i, j in group:
                        covered.update(range(i, j))
                if len(covered) < n:
                    remainder = c._intern(tuple(term for i, term in enumerate(terms)
                                                if i not in covered))
                    pieces.append(visit(remainder, budget - 1))
                candidate = pieces[0] if len(pieces) == 1 else self._node("sum", tuple(pieces))
                consider(candidate)

            def partition(cuts):
                groups = defaultdict(list)
                for i, j in pairwise((0, *cuts, n)):
                    if 2 <= j - i <= n // 2:
                        groups[signature(i, j)].append((i, j))
                build(groups.values())

            if n >= 4:
                # Break at the largest actual gaps, independent of exponent depth.
                for threshold in sorted(set(gaps), reverse=True)[:5]:
                    partition(tuple(i + 1 for i, gap in enumerate(gaps)
                                    if c.compare(gap, threshold) >= 0))

                # Prefix-function borders find any complete periodic token string.
                tokens = tuple((gap, left * right)
                               for gap, left, right in zip(gaps, digits, digits[1:]))
                prefix = [0] * len(tokens)
                for i in range(1, len(tokens)):
                    j = prefix[i - 1]
                    while j and tokens[i] != tokens[j]:
                        j = prefix[j - 1]
                    prefix[i] = j + (tokens[i] == tokens[j])
                border = prefix[-1]
                while border:
                    width = len(tokens) - border
                    if width >= 2 and n % width == 0:
                        partition(tuple(range(width, n, width)))
                    border = prefix[border - 1]

                # Allow motifs in the middle of unrelated terms. For a fixed
                # width, earliest-finish selection maximizes disjoint occurrences.
                for width in (2, 3, 4, 6, 8, 12, 16):
                    if n < 2 * width:
                        break
                    groups = defaultdict(list)
                    for i in range(n - width + 1):
                        group = groups[signature(i, i + width)]
                        if not group or group[-1][1] <= i:
                            group.append((i, i + width))
                    candidates = sorted((g for g in groups.values() if len(g) > 1),
                                        key=len, reverse=True)
                    for group in candidates[:2]:
                        build((group,))

            self._node(best.kind, best.children, best.payload,
                       small=raw._small, wide=raw._small is None)
            return best

        return visit(marking, rounds)

    def _canonicalize(self, value, remaining):
        c = self.circuit

        def spend(cost):
            if cost > remaining[0]:
                raise _BudgetExceeded
            remaining[0] -= cost

        def product(a, b):
            spend(len(a.terms) * len(b.terms))
            return multiply_fast(c, a, b)

        @cache
        def convert(expr):
            spend(1)
            kind = expr.kind
            if kind == "naf":
                terms = tuple((convert(child), -1 if expr.payload >> i & 1 else 1)
                              for i, child in enumerate(expr.children))
                spend(len(terms))
                return c._intern(terms)
            if kind == "neg":
                return c.subtract(c.zero, convert(expr.children[0]))
            if kind == "sum":
                result = c.zero
                for child in expr.children:
                    part = convert(child)
                    spend(len(result.terms) + len(part.terms))
                    result = c.add(result, part)
                return result
            if kind == "shift":
                body, exponent = map(convert, expr.children)
                spend(len(body.terms))
                return c._intern(tuple((c.add(e, exponent), d) for e, d in body.terms))
            if kind == "geo":
                spend(expr.payload)
                stride = convert(expr.children[0])
                exponent = c.zero
                result = []
                for _ in range(expr.payload):
                    result.append((exponent, 1))
                    exponent = c.add(exponent, stride)
                return c._intern(tuple(result))
            if kind == "mul":
                return product(*map(convert, expr.children))
            if kind == "repeat":
                base = convert(expr.children[0])
                count = expr.payload
                result = c.one
                while count:
                    if count & 1:
                        result = product(result, base)
                    count >>= 1
                    if count:
                        base = product(base, base)
                return result
            raise _BudgetExceeded

        return convert(value)

    def canonicalize(self, value: Node, *, max_work: int = 20_000) -> Signed:
        """Obtain a pure marking within a term-work budget, or raise ValueError."""
        self._check(value)
        try:
            return self._canonicalize(value, [max_work])
        except _BudgetExceeded as error:
            raise ValueError("pure signed normalization exceeded its budget or needs an odd quotient") from error

    def _geometric_product(self, a, b):
        if a.kind != "geo" or b.kind != "geo":
            return None
        c = self.circuit
        try:
            left = self.canonicalize(a.children[0], max_work=3000)
            right = self.canonicalize(b.children[0], max_work=3000)
        except ValueError:
            return None
        if multiply_fast(c, left, c.integer(a.payload)) is right:
            return self._node("geo", a.children, a.payload * b.payload)
        if multiply_fast(c, right, c.integer(b.payload)) is left:
            return self._node("geo", b.children, a.payload * b.payload)
        return None

    def multiply(self, a: Node, b: Node, *, max_pairs: int = 0) -> Node:
        self._check(a, b)
        if a is self.zero or b is self.zero:
            return self.zero
        for quotient, factor in ((a, b), (b, a)):
            if quotient.kind == "quotient":
                try:
                    raw_factor = self.canonicalize(factor, max_work=2000)
                except ValueError:
                    continue
                if raw_factor is self.circuit.integer(quotient.payload):
                    return quotient.children[0]
        a_base, a_count = (a.children[0], a.payload) if a.kind == "repeat" else (a, 1)
        b_base, b_count = (b.children[0], b.payload) if b.kind == "repeat" else (b, 1)
        if a_base is b_base and a_base is not self.one:
            product = self._repeat(a_base, a_count + b_count)
        else:
            product = self._geometric_product(a, b) or self._join(a, b)
        if max_pairs <= 0:
            return product
        try:
            normalized = self._canonicalize(product, [max_pairs])
        except _BudgetExceeded:
            return product
        compact = self.pack(normalized)
        return compact if len(encode(compact)) < len(encode(product)) else product

    def _remove_factor(self, value, factor):
        if value is factor:
            return self.one
        if not value._filter & factor._marker:
            return None
        if value.kind == "repeat":
            base = value.children[0]
            remainder = self._remove_factor(base, factor)
            if remainder is not None:
                return self._join(remainder, self._repeat(base, value.payload - 1))
        if value.kind == "mul":
            left, right = value.children
            remainder = self._remove_factor(right, factor)
            if remainder is not None:
                return self._join(left, remainder)
            remainder = self._remove_factor(left, factor)
            if remainder is not None:
                return self._join(remainder, right)
        return None

    def _shift_parts(self, value):
        if value.kind == "shift":
            return value.children
        if value.kind == "naf" and len(value.children) == 1 and not value.payload:
            return self.one, value.children[0]
        return value, self.zero

    def divide(self, value: Node, divisor: Node | int, *, max_terms: int = 256,
               periods=None) -> Node:
        """Cancel proven factors, recognize structured quotients or check a small divisor."""
        c = self.circuit
        self._check(value)
        if not isinstance(max_terms, int) or max_terms < 0:
            raise ValueError("max_terms must be a nonnegative integer")
        if isinstance(divisor, int):
            if divisor == 0:
                raise ValueError("division by zero")
            explicit = divisor
            divisor = self.view(c.integer(divisor))
        else:
            self._check(divisor)
            explicit = divisor._small
        if divisor is self.zero:
            raise ValueError("division by zero")
        if value is self.zero:
            return self.zero
        if divisor is self.one:
            return value
        direct = self._remove_factor(value, divisor)
        if direct is not None:
            return direct

        a, a_shift = self._shift_parts(value)
        b, b_shift = self._shift_parts(divisor)
        if a is b or b is self.one:
            try:
                a_position = self.canonicalize(a_shift, max_work=6000)
                b_position = self.canonicalize(b_shift, max_work=6000)
            except ValueError:
                a_position = b_position = None
            if a_position is not None:
                gap = c.subtract(a_position, b_position)
                if c.compare(gap, c.zero) >= 0:
                    return self._shift(self.one if a is b else a, self.pack(gap))

        if value.kind == "geo" and divisor.kind == "geo":
            if value.children == divisor.children and value.payload % divisor.payload == 0:
                stride = self.canonicalize(value.children[0], max_work=6000)
                step = multiply_fast(c, stride, c.integer(divisor.payload))
                count = value.payload // divisor.payload
                return self._node("geo", (self.pack(step),), count)

        # Match independently encoded signed operands up to a shift and sign.
        numerator_raw = divisor_raw = None
        if max_terms and (explicit is None or value.kind != "mul"):
            try:
                numerator_raw = self.canonicalize(value, max_work=12_000)
                divisor_raw = self.canonicalize(divisor, max_work=12_000)
            except ValueError:
                numerator_raw = divisor_raw = None
        if numerator_raw is not None:
            nterms, dterms = numerator_raw.terms, divisor_raw.terms
            if not dterms:
                raise ValueError("division by zero")
            if not nterms:
                return self.zero
            gap = c.subtract(nterms[0][0], dterms[0][0])
            if c.compare(gap, c.zero) >= 0:
                sign = nterms[0][1] * dterms[0][1]
                if len(nterms) == len(dterms) and all(
                    nd == sign * dd and c.subtract(ne, de) is gap
                    for (ne, nd), (de, dd) in zip(nterms, dterms)
                ):
                    quotient = self._shift(self.one, self.pack(gap))
                    return quotient if sign == 1 else self._neg(quotient)
                if len(dterms) == 1:
                    shifted = self.pack(c.divide_power(numerator_raw, dterms[0][0]))
                    return shifted if dterms[0][1] > 0 else self._neg(shifted)
            if (explicit in (3, -3) and len(nterms) == 2
                    and nterms[0][0] is c.zero and nterms[0][1] == -nterms[1][1]):
                if nterms[1][0].terms[0][0] is c.zero:
                    raise ValueError("division is not exact")
                quotient = self._node("quotient", (value,), 3)
                return quotient if explicit > 0 else self._neg(quotient)
            try:
                return self.pack(divide_fast(c, numerator_raw, divisor_raw, max_terms=max_terms))
            except QuotientWorkLimit:
                pass

        if explicit is None:
            try:
                bound = divisor_raw or self.canonicalize(divisor, max_work=6000)
                explicit = _bounded_integer(bound, 512)
            except ValueError as error:
                raise ValueError("division requires a recoverable factor or an explicit divisor") from error
        if explicit == 0:
            raise ValueError("division by zero")
        magnitude = abs(explicit)
        if residue(value, magnitude, periods=periods):
            raise ValueError("division is not exact")
        if value.kind == "quotient":
            quotient = self._node("quotient", value.children, value.payload * magnitude)
        else:
            quotient = self._node("quotient", (value,), magnitude)
        return quotient if explicit > 0 else self._neg(quotient)

    def stats(self, value):
        self._check(value)
        return stats(value)

    def encode(self, value: Node) -> bytes:
        self._check(value)
        return encode(value)

    def decode(self, data: bytes) -> Node:
        """Load an encoding from this kernel, including its arithmetic facts."""
        nodes = []
        for kind, indices, payload, wide in decode_records(data):
            children = tuple(nodes[i] for i in indices)
            small = _known_small(kind, children, payload)
            if wide and small is not None and abs(small) <= 1 << 62:
                raise ValueError("inconsistent adaptive size flag")
            nodes.append(self._node(kind, children, payload, small=small, wide=wide))
        return nodes[-1]

    def residue(self, value, modulus, *, periods=None):
        self._check(value)
        return residue(value, modulus, periods=periods)
