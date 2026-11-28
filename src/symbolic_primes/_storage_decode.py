"""Read signed-power storage rows into the canonical arithmetic circuit."""

from .symbolic import Circuit, Signed
from . import _accelerator


def decode_storage(circuit: Circuit, data: bytes) -> Signed:
    """Read one reachable H2N1 graph without materializing its values."""
    if not isinstance(data, bytes) or not data.startswith(b"H2N1"):
        raise ValueError("invalid signed circuit encoding")
    if _accelerator.native is not None and len(data) >= 256:
        row = _accelerator.native.decode_dense(data)
        if row is not None:
            return circuit._from_dense_row(row)
    offset = 4
    size = len(data)
    max_bits = (2 * size).bit_length()

    def unsigned():
        nonlocal offset
        value = shift = 0
        while offset < size:
            octet = data[offset]
            offset += 1
            value |= (octet & 127) << shift
            if octet < 128:
                if shift and not octet:
                    raise ValueError("overlong signed circuit integer")
                return value
            shift += 7
            if shift >= max_bits:
                raise ValueError("signed circuit integer exceeds format bounds")
        raise ValueError("truncated signed circuit encoding")

    count = unsigned()
    if not 0 < count <= size - offset:
        raise ValueError("invalid signed circuit node count")
    nodes = []
    defined = set()
    referenced = bytearray(count)
    for row in range(count):
        length = unsigned()
        if length > size - offset - (count - row - 1):
            raise ValueError("invalid signed circuit node arity")
        coefficients = {}
        previous = None
        for _ in range(length):
            edge = unsigned()
            child = edge >> 1
            if child >= row:
                raise ValueError("signed circuit reference must point to an earlier node")
            exponent = nodes[child]
            if exponent.terms and exponent.terms[-1][1] < 0:
                raise ValueError("negative exponent in signed circuit")
            if previous is not None and circuit.compare(previous, exponent) >= 0:
                raise ValueError("signed circuit exponents must be strictly increasing")
            coefficients[exponent] = -1 if edge & 1 else 1
            referenced[child] = 1
            previous = exponent

        value = circuit._from_coefficients(coefficients, False)
        if value in defined:
            raise ValueError("duplicate value in signed circuit")
        defined.add(value)
        nodes.append(value)

    if offset != size:
        raise ValueError("trailing data after signed circuit")
    # Every earlier node must have a path through later parents to the root.
    if not all(referenced[:-1]):
        raise ValueError("unreachable node in signed circuit")
    return nodes[-1]
