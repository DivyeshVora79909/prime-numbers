"""A whole-graph byte lower bound for wide roots at ordinary bit positions."""

from math import inf


def _ordinary_highest(value):
    # Bound the column DP by both absolute span and density. Symbolic
    # tower positions and very sparse wide graphs use the gap-skipping search.
    if len(value.terms) >= 16:
        highest = value.terms[-1][0]._small
        if highest is not None and highest <= min(8192, max(511, 16 * len(value.terms))):
            return highest
    return None


def _binary_columns(value, horizon):
    """Read the magnitude as bit columns without materializing the integer."""
    columns = [0] * (horizon + 1)
    sign = value.terms[-1][1]
    for exponent, digit in value.terms:
        columns[exponent._small] = sign * digit
    carry = 0
    for i, digit in enumerate(columns):
        coefficient = digit + carry
        columns[i] = coefficient % 2
        carry = coefficient // 2
    return columns


def _naf_weights(limit):
    weights = [0, 1]
    for x in range(2, limit + 1):
        half = x // 2
        weights.append(1 + min(weights[half], weights[half + 1]) if x % 2 else weights[half])
    return weights[: limit + 1]


def _restricted_weights(limit, cutoff, allowed, tail):
    """Minimum signed weights with specified allowed positions through cutoff."""
    above = tail
    for position in range(cutoff, -1, -1):
        scale = 1 << position
        largest = (limit + scale - 1) // scale
        row = []
        available = allowed & scale
        for x in range(largest + 1):
            half = x // 2
            if x % 2:
                row.append(1 + min(above[half], above[half + 1]) if available else inf)
            else:
                row.append(above[half])
        above = row
    return above


def _root_cost(bits, costs, *, force=-1, prefer_negative=None):
    """Cheapest signed expansion; an optional required digit position."""
    no_carry, carry = 0, inf
    for e in range(len(bits) - 1, -1, -1):
        if prefer_negative is not None:
            prefer_negative[e] = carry < no_carry
        odd = costs[e] + min(no_carry, carry)
        if bits[e]:
            no_carry, carry = odd, inf if e == force else carry
        else:
            no_carry, carry = inf if e == force else no_carry, odd
    return no_carry


def _root_digits(bits, costs):
    """Recover one minimum-cost row, preferring a positive digit on ties."""
    prefer_negative = bytearray(len(bits))
    if _root_cost(bits, costs, prefer_negative=prefer_negative) == inf:
        return None
    digits = []
    carry = 0
    for e, bit in enumerate(bits):
        coefficient = bit + carry
        digit = 0
        if coefficient % 2:
            digit = -1 if prefer_negative[e] else 1
            digits.append((e, digit))
        carry = (coefficient - digit) // 2
    return tuple(digits)


def weighted_root_row(circuit, value, step):
    """Construct an exact candidate by accounting for shared low exponents."""
    highest = _ordinary_highest(value)
    if highest is None:
        return None
    step()  # One weighted carry DP over these bounded columns.
    limit = highest + 1
    cutoff = highest.bit_length() - 1
    bits = _binary_columns(value, limit)
    weights = _naf_weights(limit)
    new_row_cost = 3 if len(value.terms) >= 65 else 2
    costs = [1 if e <= cutoff else new_row_cost + w for e, w in enumerate(weights)]
    sign = value.terms[-1][1]
    return tuple((circuit.integer(e), d * sign) for e, d in _root_digits(bits, costs))


def _secondary_cost(bits, costs, secondary):
    # A unit of the primary cost outweighs every possible secondary digit.
    scale = len(bits) + 1
    encoded = [weight * scale + extra for weight, extra in zip(costs, secondary)]
    score = _root_cost(bits, encoded)
    return (inf, 0) if score == inf else divmod(score, scale)


def _first_missing(required, allowed, start, stop):
    """Partition all proper completions by their first further missing node."""
    for e in range(start, stop):
        bit = 1 << e
        yield required, allowed ^ bit, e + 1
        required |= bit


def _late_reference(bits, cutoff, ordinary, costs, cheapest, step):
    """Can DFS order rule out equality in the full low-node relaxation?"""
    if _root_cost(bits, costs, force=cutoff) == cheapest:
        return False

    limit = len(bits) - 1
    lower_positions = (1 << cutoff) - 1
    omit = _restricted_weights(limit, cutoff, lower_positions, ordinary)
    needs_reference = [int(e == cutoff or omit[e] > w) for e, w in enumerate(ordinary)]
    if _secondary_cost(bits, costs, needs_reference)[1] < 2:
        return False

    # The first magnitude that can use a digit >= cutoff in a minimum row.
    # A nonempty signed row with distinct positions cannot sum to zero.
    nonzero_tail = [inf, *ordinary[1:]]
    uses_high = _restricted_weights(limit, cutoff - 1, lower_positions, nonzero_tail)
    first = next((e for e in range(limit + 1) if uses_high[e] == ordinary[e]), None)
    if first is None or first <= cutoff:
        return False

    # Before DFS reaches a root child >= first, its children and all their
    # descendants must use positions below cutoff. Bound a family of prefix
    # sets by allowing every undecided node and charging only required nodes
    # and the distinct undecided low nodes actually selected by the root.
    no_high_tail = [0] + [inf] * limit
    count_prior = [int(cutoff < e < first) for e in range(limit + 1)]

    def all_prefixes_late(required, allowed, start):
        step()
        within = _restricted_weights(first - 1, cutoff - 1, allowed, no_high_tail)
        eligible = costs.copy()
        optimistic_count = count_prior.copy()
        eligible[cutoff] = inf  # Equality already excluded a root digit here.
        for e in range(first):
            if within[e] > ordinary[e] or (e < cutoff and not allowed & (1 << e)):
                eligible[e] = inf
            elif e < cutoff and not required & (1 << e):
                optimistic_count[e] = 1
        primary, prior = _secondary_cost(bits, eligible, optimistic_count)
        if primary > cheapest or required.bit_count() + prior >= 64:
            return True

        # Check the completion with all currently allowed nodes present.
        primary, prior = _secondary_cost(bits, eligible, count_prior)
        if primary == cheapest and allowed.bit_count() + prior < 64:
            return False
        return all(all_prefixes_late(*group) for group in _first_missing(
            required, allowed, start, cutoff
        ))

    # An index >=64 with two incoming references pays one byte more than
    # the baseline, which allowed exactly one wide reference per such index.
    return all_prefixes_late(0, lower_positions, 0)


def ordinary_storage_bound(value, ceiling, step):
    """Bound H2N1 bytes from below, given a valid incumbent of size ceiling.

    Return None outside this bounded DP's range. Each configuration charges
    a search step. The ordinary byte bound remains valid with arbitrary
    signed descendants and every acyclic DFS order allowed by the encoder.
    """
    highest = _ordinary_highest(value)
    if highest is None:
        return None
    step()
    # Past the highest bit, a pending carry must emit one digit immediately.
    # Ending it there costs strictly less than continuing, since each
    # subsequent relaxed digit has positive cost. This bounds the relaxed
    # minimum and its ties, even if the actual optimum uses higher exponents.
    limit = highest + 1
    bits = _binary_columns(value, limit)
    ordinary = _naf_weights(limit)
    cutoff = highest.bit_length() - 1
    full = (1 << (cutoff + 1)) - 1
    minimum_terms = len(value.terms)
    wide = minimum_terms >= 65
    # Minimum root arity and reachable node count (all root children plus
    # the root) determine how many bytes those two header fields require.
    header = 4 + (minimum_terms.bit_length() + 6) // 7
    header += ((minimum_terms + 1).bit_length() + 6) // 7
    if wide:
        header -= 64
    row_header = 2 if wide else 1
    # Each selected high position adds its root edge, its own arity and
    # its defining edges.
    # With >=65 root digits there must also be an incoming multi-byte edge
    # for every non-root index >=64, hence +1 per distinct node minus 64.
    costs = [1 if e <= cutoff else 1 + row_header + w for e, w in enumerate(ordinary)]
    cheapest = _root_cost(bits, costs)
    full_bound = header + sum(row_header + ordinary[t] for t in range(cutoff + 1)) + cheapest
    # Only bounds below the full-set result (or its one-byte DFS correction)
    # can affect the answer. A family proven at least this expensive is done.
    other = full_bound + int(wide and full_bound < ceiling)

    def visit_subsets(required, allowed, start):
        nonlocal other
        for present, possible, next_bit in _first_missing(required, allowed, start, cutoff + 1):
            step()
            # Higher bits remain available, so the ordinary suffix table
            # already solves all positions above the latest exclusion.
            weights = _restricted_weights(limit, next_bit - 1, possible, ordinary)
            base = header + sum(
                row_header + weights[t] for t in range(cutoff + 1) if present & (1 << t)
            )
            if base >= other:
                continue
            weights_for_root = [1 + row_header + w for w in weights]
            for t in range(cutoff + 1):
                if present & (1 << t):
                    weights_for_root[t] = 1
                elif not possible & (1 << t):
                    weights_for_root[t] = inf
            # Undecided low nodes may serve every child row for free. If
            # the root selects one, its definition must be paid for once.
            optimistic = base + _root_cost(bits, weights_for_root)
            if optimistic >= other:
                continue
            if present == possible:
                other = optimistic
                continue
            for t in range(cutoff + 1):
                if possible & (1 << t) and not present & (1 << t):
                    base += row_header + weights[t]
                    weights_for_root[t] = 1
            other = min(other, base + _root_cost(bits, weights_for_root))
            if optimistic < other:
                visit_subsets(present, possible, next_bit)

    # Every reachable nonempty DAG contains the zero node.
    visit_subsets(1, full, 1)
    if wide and full_bound < min(other, ceiling) and _late_reference(
        bits, cutoff, ordinary, costs, cheapest, step
    ):
        full_bound += 1
    return min(other, full_bound)
