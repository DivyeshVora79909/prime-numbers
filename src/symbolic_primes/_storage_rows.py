"""Enumerate equivalent signed rows by carrying over symbolic positions."""


def alternative_rows(circuit, value, max_terms, step, *, last_position=None):
    """Yield every distinct-position row within the term cap, canonical first."""
    terms = value.terms
    stack = [(0, circuit.zero, 0, None, 0)]
    while stack:
        index, position, carry, path, count = stack.pop()
        step()
        # One pending carry can remove at most one term from the NAF suffix.
        if count + len(terms) - index - bool(carry) > max_terms:
            continue

        if carry:
            if (last_position is not None and index == len(terms)
                    and circuit.compare(position, last_position) >= 0):
                continue
            position = circuit.add(position, circuit.one)
        else:
            if index == len(terms):
                row = []
                while path is not None:
                    term, path = path
                    row.append(term)
                yield tuple(reversed(row))
                continue
            position = terms[index][0]

        coefficient = carry
        if index < len(terms) and terms[index][0] is position:
            coefficient += terms[index][1]
            index += 1

        if coefficient % 2 == 0:
            carry = coefficient // 2
            stack.append((index, position, carry, path, count))
        elif count < max_terms:
            # An odd coefficient is ±1. Match it first to stop the carry.
            opposite = ((position, -coefficient), path)
            matching = ((position, coefficient), path)
            stack.append((index, position, coefficient, opposite, count + 1))
            stack.append((index, position, 0, matching, count + 1))
