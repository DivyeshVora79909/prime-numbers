"""Move frequently used exponent definitions into the first 64 wire IDs."""

from collections import Counter
from math import inf

from ._storage_bound import _root_digits


def _first_nodes(root, rows, changed=None, replacement=()):
    # Keep the root iterator lazy. Only a short prefix of a wide root is read.
    first = set()
    stack = [(root, iter(rows[root]))]
    while stack and len(first) < 64:
        node, children = stack[-1]
        edge = next(children, None)
        if edge is None:
            first.add(node)
            stack.pop()
        else:
            child = edge[0]
            if child not in first:
                row = replacement if child is changed else rows[child]
                stack.append((child, iter(row)))
    return first


def _rows_using(value, position, available):
    # Fix the sign of this reference. The remaining signed digits can only
    # use other definitions already in the graph, once per position.
    for sign in (1, -1):
        residual = value - sign * (1 << position)
        magnitude = abs(residual)
        width = magnitude.bit_length() + 1
        bits = [(magnitude >> e) & 1 for e in range(width)]
        costs = [1 if e != position and e in available else inf for e in range(width)]
        digits = _root_digits(bits, costs)
        if digits is None:
            continue
        orientation = -1 if residual < 0 else 1
        terms = [(e, d * orientation) for e, d in digits]
        terms.append((position, sign))
        yield tuple((available[e], d) for e, d in sorted(terms))


def early_reference_change(root, rows, step):
    """Suggest one exact change to the canonical descendants of a signed root.

    For graphs with at most 8,193 rows, each edge costs two bytes minus one
    when its child is in the first 64 IDs. A changed arity or changed first
    64 set can therefore be priced without serializing the other rows. For
    wider IDs this scores the first size threshold; the caller prices the
    actual candidate before adopting it.
    """
    if len(rows) <= 65 or not root.terms:
        return None
    highest = rows[root][-1][0]._small
    if highest is None:
        return None
    max_coordinate = highest.bit_length()
    available = {
        node._small: node for node in rows
        if node._small is not None and 0 <= node._small <= max_coordinate
    }
    counts = Counter(child for row in rows.values() for child, _ in row)
    first = _first_nodes(root, rows)
    late = [e for e, node in sorted(available.items()) if node not in first and counts[node] > 1]
    if not late:
        return None
    original_credit = sum(counts[node] for node in first)
    best_delta = 0
    change = None
    for parent, _ in rows[root]:
        if parent not in first:
            break
        if parent._small <= max_coordinate:
            continue
        old = rows[parent]
        old_refs = {child for child, _ in old}
        for position in late:
            for row in _rows_using(parent._small, position, available):
                step()
                refs = {child for child, _ in row}
                removed = old_refs - refs
                if any(counts[child] == 1 for child in removed):
                    continue
                # Keeping another reference to every removed child preserves
                # reachability. All new children are numerically below parent.
                after = _first_nodes(root, rows, parent, row)
                credit = sum(counts[node] for node in after)
                credit += len((refs - old_refs) & after) - len(removed & after)
                # Both rows store an ordinary exponent with <=64 available
                # positions, so their arities still take exactly one byte.
                delta = 2 * (len(row) - len(old)) + original_credit - credit
                if delta < best_delta:
                    best_delta = delta
                    change = parent, row
    return change
