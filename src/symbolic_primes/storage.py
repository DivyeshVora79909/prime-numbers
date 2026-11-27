"""Whole-graph optimization of hereditary signed-digit storage."""

from dataclasses import dataclass
from itertools import islice
from typing import Iterator

from ._storage_bound import ordinary_storage_bound, weighted_root_row
from ._storage_decode import decode_storage
from ._storage_layout import early_reference_change
from ._storage_rows import alternative_rows
from ._storage_small import prove_small_storage
from .similarity_kernel import multiply_fast
from .symbolic import Circuit, Signed, _encode_graph_rows, _encode_rows


__all__ = ("StorageResult", "decode_storage", "multiply_storage", "optimize_storage")
_Row = tuple[tuple[Signed, int], ...]
_DEFAULT_STATES = 2_048


@dataclass(frozen=True, slots=True)
class StorageResult:
    """An exact encoding and the status of its size proof."""

    data: bytes
    optimal: bool
    states: int
    lower_bound_bytes: int = 0

    @property
    def max_extra_bytes(self) -> int:
        """The most this encoding could exceed the format's minimum by."""
        return 0 if self.optimal else len(self.data) - self.lower_bound_bytes


class _WorkLimit(Exception):
    pass


@dataclass(slots=True)
class _Frame:
    value: Signed
    rows: dict[Signed, _Row]
    pending: set[Signed]
    other_row_bytes: int
    alternatives: Iterator[_Row]


def _varint_size(value: int) -> int:
    return max(1, (value.bit_length() + 6) // 7)


def _row_minimum(count: int) -> int:
    return _varint_size(count) + count


def _minimum_size(nodes: int, row_bytes: int) -> int:
    # All reachable nodes except the last (the root) have an incoming edge.
    # Counting one reference per node gives a lower bound on wide IDs.
    wide_references = 0
    first_wide_id = 64
    while first_wide_id < nodes - 1:
        wide_references += nodes - 1 - first_wide_id
        first_wide_id *= 128
    return 4 + _varint_size(nodes) + row_bytes + wide_references


def _encode_graph(root: Signed, rows: dict[Signed, _Row]) -> bytes:
    return _encode_rows(_encode_graph_rows(root, rows.__getitem__))


def _canonical_rows(root: Signed, root_row: _Row) -> dict[Signed, _Row]:
    rows = {root: root_row}
    remaining = [child for child, _ in root_row]
    while remaining:
        node = remaining.pop()
        if node not in rows:
            rows[node] = node.terms
            remaining.extend(child for child, _ in node.terms)
    return rows


def _required_rows(values, rows, pending):
    required = set()
    remaining = list(values)
    while remaining:
        value = remaining.pop()
        if value not in required and value not in rows and value not in pending:
            required.add(value)
            # The lowest occupied position is fixed by divisibility by two.
            if value.terms:
                remaining.append(value.terms[0][0])
    return required


def _seed(circuit, root, rows, best, step):
    """Try short alternative rows on nodes already paid for by the graph."""
    for value in sorted(rows, reverse=True):
        current = rows[value]
        alternatives = alternative_rows(circuit, value, len(current) + 1, step)
        for row in islice(alternatives, 16):
            if row == current or any(child is value or child not in rows for child, _ in row):
                continue
            rows[value] = row
            candidate = _encode_graph(root, rows)
            if len(candidate) < len(best):
                current, best = row, candidate
                yield best
            else:
                rows[value] = current


def optimize_storage(
    circuit: Circuit, value: Signed, *, max_states: int | None = _DEFAULT_STATES,
    target_gap_percent: int | None = None,
) -> StorageResult:
    """Minimize H2N1 size over acyclic, shared signed-digit representations.

    Each reachable integer has one definition, with at most one signed digit
    per exponent and exponents in ascending order. Row IDs follow the same
    depth-first postorder as Circuit.encode. A bounded search aims for a
    certified gap of at most 1% by default. An unlimited search defaults
    to an exact proof. The returned bound also applies on budget exhaustion.
    """
    circuit._check(value)
    if max_states is not None and (not isinstance(max_states, int) or max_states < 0):
        raise ValueError("max_states must be a nonnegative integer or None")
    if target_gap_percent is None:
        target_gap_percent = 0 if max_states is None else 1
    elif not isinstance(target_gap_percent, int) or target_gap_percent < 0:
        raise ValueError("target_gap_percent must be a nonnegative integer or None")

    seed_rows = _canonical_rows(value, value.terms)
    best = _encode_graph(value, seed_rows)
    states = 0
    initial_pending = _required_rows((value,), {}, set())
    initial_bytes = sum(_row_minimum(len(node.terms)) for node in initial_pending)
    # Along the compulsory chain, every other child of a positive row
    # lies strictly between that row and its lowest exponent. Those
    # intervals are disjoint, also for a signed root above its lowest digit.
    # Each additional child needs a distinct nonzero row and at least one edge.
    other_children = sum(len(node.terms) - 1 for node in initial_pending if node.terms)
    lower_bound = _minimum_size(len(initial_pending) + other_children,
                                initial_bytes + 2 * other_children)

    def result(*, proved=False):
        bound = len(best) if proved else lower_bound
        assert bound <= len(best)
        return StorageResult(best, bound == len(best), states, bound)

    def within_target():
        return 100 * (len(best) - lower_bound) <= target_gap_percent * lower_bound

    if lower_bound == len(best) or max_states == 0:
        return result()

    def step():
        nonlocal states
        if max_states is not None and states >= max_states:
            raise _WorkLimit
        states += 1

    seed_limit = 5_000 if max_states is None else min(5_000, max_states // 3)
    seed_states = 0

    def seed_step():
        nonlocal seed_states
        if seed_states >= seed_limit:
            raise _WorkLimit
        step()
        seed_states += 1

    layouts = [seed_rows]
    best_rows = seed_rows
    try:
        root_row = weighted_root_row(circuit, value, seed_step)
        if root_row is not None and root_row != value.terms:
            alternative = _canonical_rows(value, root_row)
            candidate = _encode_graph(value, alternative)
            layouts.insert(0, alternative)
            if len(candidate) < len(best):
                best_rows, best = alternative, candidate
    except _WorkLimit:
        pass

    # A near-optimal run reserves most of its work for finding better bytes.
    # An interrupted subset proof contributes no provisional lower bound.
    bound_start = states
    proof_quota = None
    if max_states is not None and target_gap_percent:
        proof_quota = min(64, max_states // 4)

    def bound_step():
        if proof_quota is not None and states - bound_start >= proof_quota:
            raise _WorkLimit
        step()

    try:
        bound = ordinary_storage_bound(value, len(best), bound_step)
        if bound is not None:
            lower_bound = max(lower_bound, bound)
        if lower_bound == len(best):
            return result()
    except _WorkLimit:
        if states == max_states:
            return result()

    # Reordering early definitions can save bytes even if a weighted root
    # initially cost more than the canonical one. Price both whole graphs.
    for rows in layouts:
        if seed_states >= seed_limit:
            break
        try:
            change = early_reference_change(value, rows, seed_step)
        except _WorkLimit:
            break
        if change is not None:
            node, replacement = change
            original = rows[node]
            rows[node] = replacement
            candidate = _encode_graph(value, rows)
            if len(candidate) < len(best):
                best_rows, best = rows, candidate
            else:
                rows[node] = original
        if lower_bound == len(best):
            return result()
    if within_target():
        return result()

    def accept(data):
        nonlocal best
        if len(data) < len(best):
            best = data

    try:
        if prove_small_storage(circuit, value, len(best), step, accept):
            return result(proved=True)
    except _WorkLimit:
        return result()

    if len(best_rows) > 64 and seed_states < seed_limit:
        try:
            for candidate in _seed(circuit, value, best_rows, best, seed_step):
                best = candidate
                if within_target():
                    return result()
        except _WorkLimit:
            pass

    root_horizon = circuit.bit_length(value)

    def frame(rows, pending, row_bytes):
        node_count = len(rows) + len(pending)
        if _minimum_size(node_count, row_bytes) >= len(best):
            return None
        target = max(pending)
        other_row_bytes = row_bytes - _row_minimum(len(target.terms))
        # A prospective row costs at least one byte more than its term count.
        max_terms = len(best) - _minimum_size(node_count, other_row_bytes) - 2
        options = alternative_rows(
            circuit, target, max_terms, step,
            last_position=root_horizon if target is value else None,
        )
        return _Frame(target, rows, pending - {target}, other_row_bytes, options)

    try:
        initial = frame({}, initial_pending, initial_bytes)
        stack = [] if initial is None else [initial]
        while stack:
            current = stack[-1]
            row = next(current.alternatives, None)
            if row is None:
                stack.pop()
                continue

            row_bytes = current.other_row_bytes + _row_minimum(len(row))
            known_count = len(current.rows) + len(current.pending) + 1
            if _minimum_size(known_count, row_bytes) >= len(best):
                continue
            required = _required_rows(
                (child for child, _ in row), current.rows, current.pending
            )
            if current.value in required:
                # A positive row x using any exponent >= x must use x
                # itself: omitting position x would leave at least 2^x+1.
                # Hence every legal positive dependency strictly descends.
                continue
            row_bytes += sum(_row_minimum(len(child.terms)) for child in required)
            pending = current.pending | required
            if _minimum_size(known_count + len(required), row_bytes) >= len(best):
                continue

            definitions = {**current.rows, current.value: row}
            if pending:
                following = frame(definitions, pending, row_bytes)
                if following is not None:
                    stack.append(following)
            else:
                candidate = _encode_graph(value, definitions)
                if len(candidate) < len(best):
                    best = candidate
                    if within_target():
                        return result()
    except _WorkLimit:
        return result()
    return result(proved=True)


def multiply_storage(
    circuit: Circuit, left: Signed | bytes, right: Signed | bytes,
    *, max_states: int | None = _DEFAULT_STATES, target_gap_percent: int | None = None,
) -> StorageResult:
    """Multiply signed graphs or their H2N1 bytes, then optimize the product."""
    if isinstance(left, bytes):
        left = decode_storage(circuit, left)
    if isinstance(right, bytes):
        right = decode_storage(circuit, right)
    return optimize_storage(
        circuit, multiply_fast(circuit, left, right),
        max_states=max_states, target_gap_percent=target_gap_percent,
    )
