"""DAG serialization and modular arithmetic for adaptive symbolic integers."""

from functools import cache

from .symbolic import _totient, _varint


_OPCODES = {
    "naf": 0,
    "sum": 1,
    "neg": 2,
    "shift": 3,
    "mul": 4,
    "geo": 5,
    "quotient": 6,
    "repeat": 7,
}
_PAYLOAD_KINDS = {"geo", "quotient", "repeat"}
_FIXED_ARITY = {"neg": 1, "shift": 2, "mul": 2, "geo": 1, "quotient": 1, "repeat": 1}
_KINDS = tuple(_OPCODES)


def _serialize(root) -> tuple[bytes, int, int]:
    indices = {}
    records = []
    pending = [(root, False)]
    edges = 0
    while pending:
        node, expanded = pending.pop()
        key = id(node)
        if key in indices:
            continue
        if not expanded:
            pending.append((node, True))
            pending.extend((child, False) for child in reversed(node.children))
            continue

        row = bytearray([_OPCODES[node.kind] | (128 if getattr(node, "_wide", False) else 0)])
        _varint(row, len(node.children))
        if node.kind in _PAYLOAD_KINDS:
            _varint(row, node.payload)
        for i, child in enumerate(node.children):
            index = indices[id(child)]
            if node.kind == "naf":
                index = 2 * index + (node.payload >> i & 1)
            _varint(row, index)
        indices[key] = len(records)
        records.append(row)
        edges += len(node.children)

    data = bytearray(b"H2A1")
    _varint(data, len(records))
    data.extend(b"".join(records))
    return bytes(data), len(records), edges


def encode(expr) -> bytes:
    """Encode each distinct reachable node once, with the root last."""
    return _serialize(expr)[0]


def stats(expr) -> dict:
    data, nodes, edges = _serialize(expr)
    return {"bytes": len(data), "nodes": nodes, "edges": edges}


def decode_records(data: bytes):
    """Parse records emitted by encode, checking topology and record shapes."""
    if not isinstance(data, bytes) or data[:4] != b"H2A1":
        raise ValueError("invalid adaptive circuit encoding")
    index = 4

    def unsigned():
        nonlocal index
        result = shift = 0
        while index < len(data):
            octet = data[index]
            index += 1
            result |= (octet & 127) << shift
            if octet < 128:
                return result
            shift += 7
        raise ValueError("truncated adaptive circuit encoding")

    count = unsigned()
    if not 0 < count <= len(data) - index:
        raise ValueError("invalid adaptive node count")
    for row in range(count):
        if index >= len(data):
            raise ValueError("truncated adaptive circuit encoding")
        tag = data[index]
        index += 1
        code = tag & 127
        if code >= len(_KINDS):
            raise ValueError("unknown adaptive circuit opcode")
        kind = _KINDS[code]
        length = unsigned()
        if length > len(data) - index or length != _FIXED_ARITY.get(kind, length):
            raise ValueError("invalid adaptive node arity")
        if kind == "sum" and length < 2:
            raise ValueError("invalid adaptive sum arity")
        payload = unsigned() if kind in _PAYLOAD_KINDS else None
        if payload is not None and payload < (2 if kind == "repeat" else 1):
            raise ValueError("invalid adaptive node payload")
        if kind == "naf":
            payload = 0
        children = []
        for position in range(length):
            child = unsigned()
            if kind == "naf":
                payload |= (child & 1) << position
                child >>= 1
            if child >= row:
                raise ValueError("adaptive reference must point to an earlier node")
            children.append(child)
        yield kind, tuple(children), payload, bool(tag & 128)
    if index != len(data):
        raise ValueError("trailing data after adaptive circuit")


def _geometric_mod(ratio: int, count: int, modulus: int) -> int:
    running_power, running_sum = 1, 0
    for bit in bin(count)[2:]:
        running_sum = running_sum * (1 + running_power) % modulus
        running_power = running_power * running_power % modulus
        if bit == "1":
            running_sum = (running_sum + running_power) % modulus
            running_power = running_power * ratio % modulus
    return running_sum


def _clip(value: int, cap: int) -> int:
    return max(-cap, min(value, cap))


def _capped_power(base: int, exponent: int, cap: int) -> int:
    factor = abs(base)
    if factor <= 1:
        return pow(base, exponent)
    sign = -1 if base < 0 and exponent % 2 else 1
    if exponent >= cap.bit_length():
        return sign * cap
    magnitude = 1
    while exponent:
        if exponent & 1:
            magnitude = min(magnitude * factor, cap)
        exponent >>= 1
        factor = min(factor * factor, cap)
    return sign * magnitude


def residue(expr, modulus: int, *, periods=None) -> int:
    """Evaluate a residue without constructing the represented integer."""
    if not isinstance(modulus, int) or modulus <= 0:
        raise ValueError("modulus must be a positive int")
    periods = {} if periods is None else periods
    residues = {}
    powers = {}
    clipped = {}

    @cache
    def period(odd):
        if odd in periods:
            value = periods[odd]
            if not isinstance(value, int) or value <= 0 or pow(2, value, odd) != 1:
                raise ValueError("invalid supplied period for 2")
            return value
        if odd >= 1 << 32:
            raise ValueError("recursive exponent needs a supplied period for this odd modulus")
        return _totient(odd)

    def bounded(node, cap):
        key = id(node), cap
        if key not in clipped:
            small = getattr(node, "_small", None)
            clipped[key] = _clip(small, cap) if small is not None else bounded_node(node, cap)
        return clipped[key]

    def bounded_nonnegative(node, cap):
        if getattr(node, "_wide", False) and cap <= 1 << 62:
            return cap
        value = bounded(node, cap)
        if value < 0:
            raise ValueError("expected a nonnegative exponent")
        return value

    def bounded_sum(children, cap):
        # Recompute after adding enough headroom for the exactly known terms.
        # Opposite-signed saturated terms cannot be subtracted soundly.
        limit = cap
        small_cancellation_limit = max(cap, 1 << 16)
        while True:
            values = [bounded(child, limit) for child in children]
            if -limit in values and limit in values:
                if limit >= small_cancellation_limit:
                    raise ValueError("threshold needs a canonical view of a signed cancellation")
                limit = min(2 * limit, small_cancellation_limit)
                continue
            known = sum(value for value in values if abs(value) < limit)
            if -limit not in values and limit not in values:
                return _clip(known, cap)
            sign = -1 if -limit in values else 1
            needed = max(cap, cap - sign * known)
            if needed <= limit:
                return sign * cap
            limit = needed

    def bounded_node(node, cap):
        kind, children, payload = node.kind, node.children, node.payload
        if kind == "naf":
            if not children:
                return 0
            # Nonadjacency bounds the tail's magnitude below a third of the
            # leading term. An exponent this large establishes saturation.
            cutoff = cap.bit_length() + 1
            if bounded_nonnegative(children[-1], cutoff) == cutoff:
                return -cap if payload >> (len(children) - 1) & 1 else cap
            return _clip(sum(
                (-1 if payload >> i & 1 else 1) << bounded_nonnegative(child, cutoff)
                for i, child in enumerate(children)
            ), cap)
        if kind == "sum":
            return bounded_sum(children, cap)
        if kind == "neg":
            return -bounded(children[0], cap)
        if kind == "shift":
            value = bounded(children[0], cap)
            if not value:
                return 0
            exponent = bounded_nonnegative(children[1], cap.bit_length())
            return _clip(value << exponent, cap)
        if kind == "mul":
            value = 1
            for child in children:
                value = _clip(value * bounded(child, cap), cap)
                if not value:
                    break
            return value
        if kind == "repeat":
            return _capped_power(bounded(children[0], max(cap, 2)), payload, cap)
        if kind == "geo":
            if payload <= 1:
                return min(payload, cap)
            stride = bounded_nonnegative(children[0], cap.bit_length())
            if not stride:
                return min(payload, cap)
            if stride * (payload - 1) >= cap.bit_length():
                return cap
            return min(((1 << (stride * payload)) - 1) // ((1 << stride) - 1), cap)
        if kind == "quotient":
            return bounded(children[0], cap * payload) // payload
        raise ValueError(f"unsupported symbolic node {kind!r}")

    def power_mod(exponent, m):
        key = id(exponent), m
        if key not in powers:
            small = getattr(exponent, "_small", None)
            if small is None:
                powers[key] = recursive_power_mod(exponent, m)
            else:
                if small < 0:
                    raise ValueError("expected a nonnegative exponent")
                powers[key] = pow(2, small, m)
        return powers[key]

    def recursive_power_mod(exponent, m):
        twos = (m & -m).bit_length() - 1
        two_power = 1 << twos
        odd = m >> twos
        mod_two = 0
        if twos:
            power = bounded_nonnegative(exponent, twos)
            if power < twos:
                mod_two = 1 << power
        if odd == 1:
            return mod_two
        mod_odd = pow(2, evaluate(exponent, period(odd)), odd)
        if not twos:
            return mod_odd
        return mod_two + two_power * ((mod_odd - mod_two) * pow(two_power, -1, odd) % odd)

    def evaluate(node, m):
        if m == 1:
            return 0
        small = getattr(node, "_small", None)
        if small is not None:
            return small % m
        key = id(node), m
        if key not in residues:
            residues[key] = evaluate_node(node, m)
        return residues[key]

    def evaluate_node(node, m):
        kind, children, payload = node.kind, node.children, node.payload
        if kind == "naf":
            return sum(
                (-1 if payload >> i & 1 else 1) * power_mod(child, m)
                for i, child in enumerate(children)
            ) % m
        if kind == "sum":
            return sum(evaluate(child, m) for child in children) % m
        if kind == "neg":
            return -evaluate(children[0], m) % m
        if kind == "shift":
            body = evaluate(children[0], m)
            return body * power_mod(children[1], m) % m if body else 0
        if kind == "mul":
            value = 1
            for child in children:
                value = value * evaluate(child, m) % m
                if not value:
                    break
            return value
        if kind == "repeat":
            return pow(evaluate(children[0], m), payload, m)
        if kind == "geo":
            if payload <= 1:
                return payload % m
            return _geometric_mod(power_mod(children[0], m), payload, m)
        if kind == "quotient":
            result, remainder = divmod(evaluate(children[0], m * payload), payload)
            if remainder:
                raise ValueError("division is not exact")
            return result
        raise ValueError(f"unsupported symbolic node {kind!r}")

    return evaluate(expr, modulus)
