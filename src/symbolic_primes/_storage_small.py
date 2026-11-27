"""Exact shared-row search when every improving H2N1 graph has one-byte fields."""

from functools import cache

from .symbolic import _encode_graph_rows, _encode_rows


def prove_small_storage(circuit, root, ceiling, step, accept):
    """Report whether the one-byte proof applies; emit each improvement."""
    # A graph with n >= 66 has at least one wide reference. Even ignoring
    # that premium it costs 5+n bytes for header/rows and at least m+n-2
    # edges: the root needs m digits and every other nonzero row needs one.
    if ceiling > 135 + len(root.terms):
        return False
    best = ceiling
    # A degree-r row needs r distinct children, at least r-1 of them nonzero.
    degree_cap = min(64, (ceiling - 6) // 3)
    if len(root.terms) > degree_cap:
        return True

    @cache
    def columns(value):
        result = []
        terms = value.terms
        for i, (position, digit) in enumerate(terms):
            gap = None
            if i + 1 < len(terms):
                distance = circuit.subtract(terms[i + 1][0], position)._small
                if distance is not None and distance - 1 <= degree_cap:
                    gap = distance
            result.append((position, digit, gap))
        return tuple(result)

    @cache
    def required(value):
        # A carry traversing a gap must emit a nonzero digit at every vacant
        # position. After an unbridgeable gap the first position is mandatory.
        found = {value}
        todo = [value]
        while todo:
            previous_gap = None
            for position, _, gap in columns(todo.pop()):
                if previous_gap is None and position not in found:
                    found.add(position)
                    todo.append(position)
                previous_gap = gap
        return frozenset(found)

    @cache
    def gap_positions(start, gap):
        positions = []
        for _ in range(gap - 1):
            start = circuit.add(start, circuit.one)
            positions.append(start)
        return tuple(positions)

    def rows(value, available, allowance):
        source = columns(value)

        def edge_cost(position):
            # Each row has distinct children. New children cost at least one
            # row byte and their canonical minimum digit count.
            return 1 if position in available else 2 + len(position.terms)

        @cache
        def moves(index, carry):
            position, digit, gap = source[index]
            coefficient = digit + carry
            if coefficient == 0:
                return [(0, (), 0)]
            if abs(coefficient) == 1:
                sign = coefficient
                options = [(0, ((position, sign),), edge_cost(position))]
                crossing = ((position, -sign),)
            else:
                sign = coefficient // 2
                following = circuit.add(position, circuit.one)
                options = [(0, ((following, sign),), edge_cost(following))]
                crossing = ()

            # Stopping any later inside a vacant gap adds edges to a row
            # with identical value and a strictly smaller support. With
            # single-byte references that replacement always saves space.
            if gap is not None:
                bridge = gap_positions(position, gap)
                extended = crossing + tuple((p, -sign) for p in bridge)
                options.append((sign, extended, sum(edge_cost(p) for p, _ in extended)))
            return options

        @cache
        def suffix(index, carry):
            step()
            if index == len(source):
                return 0
            return min(cost + suffix(index + 1, outgoing)
                       for outgoing, _, cost in moves(index, carry))

        def visit(index, carry, prefix, cost):
            step()
            if cost > allowance or len(prefix) > degree_cap:
                return
            # A two/three-state carry DP bounds the as-yet-unwritten suffix.
            # For short rows direct enumeration avoids the extra setup.
            if len(source) >= 5 and cost + suffix(index, carry) > allowance:
                return
            if index == len(source):
                yield prefix
                return
            for outgoing, segment, extra in moves(index, carry):
                yield from visit(index + 1, outgoing, prefix + segment, cost + extra)

        yield from visit(0, 0, (), 0)

    definitions = {}
    cheapest_prefix = {}

    def search(pending, paid):
        nonlocal best
        step()
        lower = paid + sum(1 + len(value.terms) for value in pending)
        if lower >= best or cheapest_prefix.get(pending, best) <= paid:
            return
        cheapest_prefix[pending] = paid
        if not pending:
            data = _encode_rows(_encode_graph_rows(root, definitions.__getitem__))
            best = len(data)
            accept(data)
            return

        # Above the highest canonical position a carry must occupy *every*
        # position until it stops. Any positive row v reaching e >= v
        # would therefore reference v itself. For every acyclic graph the
        # largest pending positive row can be settled once and forgotten.
        value = root if root in pending else max(pending)
        available = pending - {value}
        allowance = best - 1 - (lower - len(value.terms))
        positive = value.terms and value.terms[-1][1] > 0
        for row in rows(value, available, allowance):
            if positive and row and circuit.compare(row[-1][0], value) >= 0:
                continue
            unseen = (required(p) for p, _ in row if p not in available)
            following = available.union(*unseen)
            definitions[value] = row
            search(following, paid + 1 + len(row))
            del definitions[value]

    search(required(root), 5)
    return True
