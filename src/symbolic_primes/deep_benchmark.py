"""Reproducible benchmarks for random recursive and prime-product operands."""

import argparse
import json
import platform
from dataclasses import dataclass
from functools import cache, cmp_to_key
from random import Random
from statistics import median

from .benchmark import measure, prime_product_mod, prime_sample
from .symbolic import Circuit, exact_divide, residue


MODULI = (97, 257, 1009, 3)
DEFAULT_CONFIGS = ((3, 2), (3, 4), (3, 8), (5, 2), (5, 3), (8, 2))
DEFAULT_SEEDS = (51, 342, 20261128)


@dataclass(frozen=True, slots=True, eq=False)
class RawMarking:
    """Independent fixture syntax with descending nonadjacent power positions."""

    terms: tuple[tuple["RawMarking | int", int], ...]


class FixtureBuilder:
    def __init__(self, width):
        self.width = width

        @cache
        def compare(left, right):
            if left is right:
                return 0
            if isinstance(left, int) and isinstance(right, int):
                return (left > right) - (left < right)
            if isinstance(left, int):
                return -1
            if isinstance(right, int):
                return 1
            for (lexp, lsign), (rexp, rsign) in zip(left.terms, right.terms):
                order = compare(lexp, rexp)
                if order < 0:
                    return -rsign
                if order > 0:
                    return lsign
                if lsign != rsign:
                    return lsign - rsign
            if len(left.terms) > len(right.terms):
                return left.terms[len(right.terms)][1]
            if len(right.terms) > len(left.terms):
                return -right.terms[len(left.terms)][1]
            return 0

        self.compare = compare
        self.key = cmp_to_key(compare)

    def marking(self, children, rng):
        children = sorted(children, key=self.key, reverse=True)
        terms = [(children[0], 1)]
        terms.extend((child, rng.choice((-1, 1))) for child in children[1:])
        return RawMarking(tuple(terms))

    def tree(self, depth, rng):
        if depth == 1:
            return self.marking(rng.sample(range(128, 516, 4), self.width), rng)
        children = []
        while len(children) < self.width:
            candidate = self.tree(depth - 1, rng)
            if all(self.compare(candidate, child) for child in children):
                children.append(candidate)
        return self.marking(children, rng)

    def template_pool(self, depth, rng):
        leaves = range(128, 516, 4)
        pool = rng.sample(leaves, max(8, self.width + 2))
        for _ in range(depth - 1):
            next_pool = []
            while len(next_pool) < len(pool):
                candidate = self.marking(rng.sample(pool, self.width), rng)
                if all(self.compare(candidate, choice) for choice in next_pool):
                    next_pool.append(candidate)
            pool = next_pool
        return pool

    def pair(self, depth, seed, shared_templates):
        first = Random(seed)
        second = Random(seed ^ 0x2BE16B)
        if shared_templates:
            pool = self.template_pool(depth, Random(seed ^ 0x7F3021))
            left = self.marking(first.sample(pool, self.width), first)
            right = self.marking(second.sample(pool, self.width), second)
        else:
            left, right = self.tree(depth, first), self.tree(depth, second)
        if first.randrange(2):
            left = RawMarking(tuple((exponent, -sign) for exponent, sign in left.terms))
        if second.randrange(2):
            right = RawMarking(tuple((exponent, -sign) for exponent, sign in right.terms))
        return left, right


@cache
def doubling_cycle(modulus):
    """Find the exact residue cycle by inspection, with no period factoring."""
    observed = {}
    path = []
    value = 1 % modulus
    while value not in observed:
        observed[value] = len(path)
        path.append(value)
        value = (2 * value) % modulus
    return path, observed[value]


def raw_residues(value, moduli):
    @cache
    def evaluate(node, modulus):
        if isinstance(node, int):
            return node % modulus
        if modulus == 1:
            return 0
        path, start = doubling_cycle(modulus)
        period = len(path) - start
        result = 0
        for exponent, sign in node.terms:
            if isinstance(exponent, int):
                result += sign * pow(2, exponent, modulus)
            else:
                # All generated markings exceed 2**127 and the cycle's
                # preperiod is at most the bit length of the small modulus.
                index = start + (evaluate(exponent, period) - start) % period
                result += sign * path[index]
        return result % modulus

    return {modulus: evaluate(value, modulus) for modulus in moduli}


def import_pair(circuit, first, second):
    @cache
    def convert(node):
        if isinstance(node, int):
            return circuit.integer(node)
        value = circuit.zero
        for exponent, sign in node.terms:
            power = circuit.power(convert(exponent))
            if sign > 0:
                value = circuit.add(value, power)
            else:
                value = circuit.subtract(value, power)
        return value

    return convert(first), convert(second)


def raw_graph_stats(first, second):
    def descendants(root):
        found = set()
        pending = [root]
        while pending:
            current = pending.pop()
            if isinstance(current, int) or current in found:
                continue
            found.add(current)
            pending.extend(exponent for exponent, _ in current.terms)
        return found

    left, right = descendants(first), descendants(second)
    return {
        "left_recursive_nodes": len(left),
        "right_recursive_nodes": len(right),
        "shared_recursive_nodes": len(left & right),
    }


def pure_deep_case(raw_left, raw_right, references, *, fast=False, canonical_divide=False):
    from .similarity_kernel import multiply_fast

    circuit = Circuit()
    (left, right), import_us = measure(
        lambda: import_pair(circuit, raw_left, raw_right)
    )
    operation = (lambda: multiply_fast(circuit, left, right)) if fast else (
        lambda: circuit.multiply(left, right)
    )
    product, multiply_us = measure(operation)
    mod_us = None
    for modulus, (left_mod, right_mod) in references.items():
        assert residue(circuit, left, modulus) == left_mod
        assert residue(circuit, right, modulus) == right_mod
        observed, elapsed = measure(lambda: residue(circuit, product, modulus))
        assert observed == left_mod * right_mod % modulus
        if modulus == 97:
            mod_us = elapsed

    divisor = circuit.integer(3)
    numerator = circuit.multiply(product, divisor)
    quotient, divide_us = measure(lambda: exact_divide(circuit, numerator, divisor))
    for modulus, (left_mod, right_mod) in references.items():
        assert residue(circuit, quotient, modulus) == left_mod * right_mod % modulus
    restored, restore_us = measure(lambda: circuit.multiply(quotient, divisor))
    assert restored is numerator

    shift = left.terms[-1][0]
    shifted = circuit.multiply(product, circuit.power(shift))
    shifted_quotient, shift_divide_us = measure(
        lambda: circuit.divide_power(shifted, shift)
    )
    assert shifted_quotient is product
    results = {
        "input_import_us": import_us,
        "left_stats": circuit.stats(left),
        "right_stats": circuit.stats(right),
        "product_stats": circuit.stats(product),
        "multiply_us": multiply_us,
        "product_mod_97_us": mod_us,
        "checked_odd_divide_us": divide_us,
        "odd_quotient_type": type(quotient).__name__,
        "multiply_back_us": restore_us,
        "shift_divide_us": shift_divide_us,
    }
    if canonical_divide:
        from .quotient_kernel import divide_fast

        # The Signed numerator carries a normalized marking without a
        # multiplication node. The kernel receives no provenance.
        recovered, canonical_divide_us = measure(
            lambda: divide_fast(circuit, product, right, max_terms=128)
        )
        assert recovered is left
        results["canonical_divide"] = {
            "divide_us": canonical_divide_us,
            "numerator_stats": circuit.stats(product),
            "denominator_stats": circuit.stats(right),
            "quotient_stats": circuit.stats(recovered),
            "exact_canonical_match": True,
        }
    return results, circuit.encode(product)


def adaptive_deep_case(raw_left, raw_right, references, *, mining, eager):
    from .adaptive import Adaptive

    circuit = Circuit()
    left, right = import_pair(circuit, raw_left, raw_right)
    adaptive = Adaptive(circuit)
    convert = adaptive.pack if mining else adaptive.view
    (left_node, right_node), convert_us = measure(
        lambda: (convert(left), convert(right))
    )
    if eager:
        operation = lambda: adaptive.multiply(left_node, right_node, max_pairs=1_000)
    else:
        operation = lambda: adaptive.multiply(left_node, right_node)
    product, multiply_us = measure(operation)
    mod_us = None
    for modulus, (left_mod, right_mod) in references.items():
        assert adaptive.residue(left_node, modulus) == left_mod
        assert adaptive.residue(right_node, modulus) == right_mod
        observed, elapsed = measure(lambda: adaptive.residue(product, modulus))
        assert observed == left_mod * right_mod % modulus
        if modulus == 97:
            mod_us = elapsed

    three = adaptive.view(circuit.integer(3))
    numerator = adaptive.multiply(product, three)
    quotient, divide_us = measure(lambda: adaptive.divide(numerator, three))
    restored, restore_us = measure(lambda: adaptive.multiply(quotient, three))
    for modulus, (left_mod, right_mod) in references.items():
        expected = left_mod * right_mod % modulus
        assert adaptive.residue(quotient, modulus) == expected
        assert adaptive.residue(restored, modulus) == 3 * expected % modulus

    return {
        "input_conversion_us": convert_us,
        "left_stats": adaptive.stats(left_node),
        "right_stats": adaptive.stats(right_node),
        "product_stats": adaptive.stats(product),
        "product_kind": product.kind,
        "multiply_us": multiply_us,
        "product_mod_97_us": mod_us,
        "structural_divide_us": divide_us,
        "quotient_kind": quotient.kind,
        "multiply_back_us": restore_us,
    }


def deep_cases(configs, seeds):
    rows = []
    for depth, width in configs:
        for shared_templates in (False, True):
            for seed in seeds:
                builder = FixtureBuilder(width)
                raw_left, raw_right = builder.pair(depth, seed, shared_templates)
                left_mod = raw_residues(raw_left, MODULI)
                right_mod = raw_residues(raw_right, MODULI)
                references = {
                    modulus: (left_mod[modulus], right_mod[modulus])
                    for modulus in MODULI
                }
                pure, encoded = pure_deep_case(raw_left, raw_right, references)
                fast, fast_encoded = pure_deep_case(
                    raw_left, raw_right, references, fast=True,
                    canonical_divide=(seed == seeds[0]),
                )
                assert fast_encoded == encoded
                rows.append(
                    {
                        "depth": depth,
                        "width": width,
                        "seed": seed,
                        "shared_templates": shared_templates,
                        "operand_signs": (raw_left.terms[0][1], raw_right.terms[0][1]),
                        "raw_input_stats": raw_graph_stats(raw_left, raw_right),
                        "pure": pure,
                        "pure_fast": fast,
                        "adaptive_without_mining": adaptive_deep_case(
                            raw_left, raw_right, references, mining=False, eager=False
                        ),
                        "adaptive_packed": adaptive_deep_case(
                            raw_left, raw_right, references, mining=True, eager=False
                        ),
                        "adaptive_eager": adaptive_deep_case(
                            raw_left, raw_right, references, mining=True, eager=True
                        ),
                    }
                )
    return rows


def motif_control(seed):
    """Place a repeated signed block at the root and one level down."""
    from .adaptive import Adaptive

    circuit = Circuit()
    adaptive = Adaptive(circuit)
    rng = Random(seed)
    positions = sorted(rng.sample(range(128, 512, 4), 6))
    signs = [rng.choice((-1, 1)) for _ in positions]
    signs[-1] = 1
    block = circuit.zero
    for position, sign in zip(positions, signs):
        term = circuit.power(circuit.integer(position))
        block = circuit.add(block, term) if sign > 0 else circuit.subtract(block, term)
    # The repeated block has fixed inner gaps and explicit separators.
    offsets = tuple(range(0, 8 * 512, 512))
    shifts = circuit.zero
    for offset in offsets:
        shifts = circuit.add(shifts, circuit.power(circuit.integer(offset)))
    repeated = circuit.multiply(block, shifts)
    nested = circuit.power(repeated)

    def reference(modulus):
        base = sum(sign * pow(2, exponent, modulus)
                   for exponent, sign in zip(positions, signs))
        return base * sum(pow(2, offset, modulus) for offset in offsets) % modulus

    results = {}
    for label, marking in (("root", repeated), ("nested", nested)):
        baseline = adaptive.view(marking)
        packed, pack_us = measure(lambda: adaptive.pack(marking))
        for modulus in MODULI:
            expected = reference(modulus)
            if label == "nested":
                cycle, start = doubling_cycle(modulus)
                period = len(cycle) - start
                expected = cycle[start + (reference(period) - start) % period]
            assert adaptive.residue(baseline, modulus) == expected
            assert adaptive.residue(packed, modulus) == expected
        results[label] = {
            "pure_stats": circuit.stats(marking),
            "unmined_stats": adaptive.stats(baseline),
            "packed_stats": adaptive.stats(packed),
            "packed_kind": packed.kind,
            "pack_us": pack_us,
        }
    return {
        "seed": seed,
        "terms_in_block": len(positions),
        "repetitions": len(offsets),
        "stride": 512,
        "results": results,
    }


def shared_prefix_control(seed):
    """Contrast pure kernels when many positions share one deep high part."""
    from .similarity_kernel import multiply_fast

    prefix_count, terms_per_prefix = 4, 32
    depth = 3
    prefix_fixture = FixtureBuilder(prefix_count).tree(depth, Random(seed))
    prefixes = tuple(exponent for exponent, _ in reversed(prefix_fixture.terms))

    def offsets(seed_for_operand):
        rng = Random(seed_for_operand)
        return tuple(
            tuple((offset, rng.choice((-1, 1)))
                  for offset in sorted(rng.sample(range(-256, 257, 4), terms_per_prefix)))
            for _ in prefixes
        )

    left_offsets, right_offsets = offsets(seed ^ 0x2A9DD), offsets(seed ^ 0xA9D2D)
    expected = {}
    for modulus in MODULI:
        power_values, start = doubling_cycle(modulus)
        period = len(power_values) - start
        prefix_residues = [raw_residues(prefix, (period,))[period] for prefix in prefixes]

        def reference(blocks):
            return sum(
                sign * power_values[start + (high + offset - start) % period]
                for high, block in zip(prefix_residues, blocks)
                for offset, sign in block
            ) % modulus

        expected[modulus] = reference(left_offsets) * reference(right_offsets) % modulus

    trials, encodings = {}, []
    for fast in (False, True):
        circuit = Circuit()
        marking, _ = import_pair(circuit, prefix_fixture, prefix_fixture)
        highs = tuple(exponent for exponent, _ in marking.terms)

        def instantiate(blocks):
            # The source positions are multiples of at least 2**128.
            # Distinct blocks remain separated after bounded offsets.
            return circuit._intern(tuple(
                (circuit.add(high, circuit.integer(offset)), sign)
                for high, block in zip(highs, blocks)
                for offset, sign in block
            ))

        left, right = instantiate(left_offsets), instantiate(right_offsets)
        operation = (lambda: multiply_fast(circuit, left, right)) if fast else (
            lambda: circuit.multiply(left, right)
        )
        product, multiply_us = measure(operation)
        for modulus, expected_value in expected.items():
            assert residue(circuit, product, modulus) == expected_value
        trials["fast" if fast else "original"] = {
            "multiply_us": multiply_us,
            "product_stats": circuit.stats(product),
        }
        encodings.append(circuit.encode(product))
    assert encodings[0] == encodings[1]
    return {
        "seed": seed,
        "input_depth": depth,
        "shared_high_parts": prefix_count,
        "terms_per_high_part": terms_per_prefix,
        "raw_term_pairs": (prefix_count * terms_per_prefix) ** 2,
        "high_part_pairs": prefix_count ** 2,
        "trials": trials,
    }


def independent_wide_control(seed, width=128):
    """Large root supports with independently generated signed exponent branches."""
    from .similarity_kernel import multiply_fast

    builder = FixtureBuilder(width=2)

    def marking(rng):
        branches = []
        while len(branches) < width:
            candidate = builder.tree(2, rng)
            if all(builder.compare(candidate, branch) for branch in branches):
                branches.append(candidate)
        ordered = sorted(branches, key=builder.key, reverse=True)
        return RawMarking(tuple((branch, rng.choice((-1, 1))) for branch in ordered))

    left_raw = marking(Random(seed ^ 0xA4389))
    right_raw = marking(Random(seed ^ 0x74656))
    left_mod = raw_residues(left_raw, MODULI)
    right_mod = raw_residues(right_raw, MODULI)
    trials = {}
    encodings = []
    for label in ("original", "fast"):
        circuit = Circuit()
        left, right = import_pair(circuit, left_raw, right_raw)
        if label == "original":
            operation = lambda: circuit.multiply(left, right)
        else:
            operation = lambda: multiply_fast(circuit, left, right)
        product, elapsed = measure(operation)
        for modulus in MODULI:
            assert residue(circuit, product, modulus) == left_mod[modulus] * right_mod[modulus] % modulus
        trials[label] = {"multiply_us": elapsed, "product_stats": circuit.stats(product)}
        encodings.append(circuit.encode(product))
    assert encodings[0] == encodings[1]
    return {"seed": seed, "input_depth": 3, "terms_per_operand": width,
            "raw_input_stats": raw_graph_stats(left_raw, right_raw), "trials": trials}


def pure_prime_cases(primes, outsider, checkpoints, *, fast):
    from .similarity_kernel import multiply_fast

    circuit = Circuit()
    product = circuit.one
    times = []
    rows = []
    for count, prime in enumerate(primes, 1):
        factor = circuit.integer(prime)
        if fast:
            operation = lambda: multiply_fast(circuit, product, factor)
        else:
            operation = lambda: circuit.multiply(product, factor)
        product, elapsed = measure(operation)
        times.append(elapsed)
        if count not in checkpoints:
            continue

        factors = primes[:count]
        expected = prime_product_mod(factors, outsider)
        observed, mod_us = measure(lambda: residue(circuit, product, outsider))
        assert observed == expected and expected != 0
        assert residue(circuit, product, prime) == 0
        removed = primes[count // 2]
        divisor = circuit.integer(removed)
        quotient, divide_us = measure(lambda: exact_divide(circuit, product, divisor))
        remainder_factors = [factor for factor in factors if factor != removed]
        deleted, deleted_query_us = measure(lambda: residue(circuit, quotient, removed))
        assert deleted == prime_product_mod(remainder_factors, removed)
        assert deleted != 0
        restored, restore_us = measure(lambda: circuit.multiply(quotient, divisor))
        assert restored is product
        highest = product.terms[-1][0]._small
        assert highest is not None
        below_leading_power = len(product.terms) > 1 and product.terms[-2][1] < 0
        bit_length = highest if below_leading_power else highest + 1
        rows.append({
            "prime_factors": count,
            "packed_binary_bytes": (bit_length + 7) // 8,
            "product_stats": circuit.stats(product),
            "interval_median_update_us": median(times),
            "nonmember_query_us": mod_us,
            "checked_odd_divide_us": divide_us,
            "query_removed_us": deleted_query_us,
            "multiply_back_us": restore_us,
            "quotient_kind": type(quotient).__name__,
        })
        times.clear()
    return rows


def adaptive_prime_cases(primes, outsider, checkpoints, *, mining):
    from .adaptive import Adaptive

    circuit = Circuit()
    adaptive = Adaptive(circuit)
    product = adaptive.one
    convert = adaptive.pack if mining else adaptive.view
    memberships = []
    update_times = []
    input_times = []
    rows = []
    for count, prime in enumerate(primes, 1):
        factor, input_us = measure(lambda: convert(circuit.integer(prime)))
        memberships.append(factor)
        product, update_us = measure(lambda: adaptive.multiply(product, factor))
        input_times.append(input_us)
        update_times.append(update_us)
        if count not in checkpoints:
            continue

        factors = primes[:count]
        expected, baseline_us = measure(lambda: prime_product_mod(factors, outsider))
        observed, mod_us = measure(lambda: adaptive.residue(product, outsider))
        assert observed == expected and expected != 0
        assert adaptive.residue(product, prime) == 0
        removed = primes[count // 2]
        divisor = memberships[count // 2]
        quotient, divide_us = measure(lambda: adaptive.divide(product, divisor))
        remainder_factors = [factor for factor in factors if factor != removed]
        deleted, deleted_query_us = measure(
            lambda: adaptive.residue(quotient, removed)
        )
        assert deleted == prime_product_mod(remainder_factors, removed)
        assert deleted != 0
        assert adaptive.residue(quotient, outsider) == prime_product_mod(
            remainder_factors, outsider
        )
        restored, restore_us = measure(lambda: adaptive.multiply(quotient, divisor))
        assert adaptive.residue(restored, outsider) == expected
        assert adaptive.residue(restored, removed) == 0
        rows.append({
            "prime_factors": count,
            "packed_expanded_bytes_upper_bound": (129 * count + 7) // 8,
            "product_stats": adaptive.stats(product),
            "product_kind": product.kind,
            "product_tree_height": product._height,
            "interval_median_input_us": median(input_times),
            "interval_median_update_us": median(update_times),
            "nonmember_query_us": mod_us,
            "explicit_factor_list_query_us": baseline_us,
            "known_factor_divide_us": divide_us,
            "query_removed_us": deleted_query_us,
            "multiply_back_us": restore_us,
            "quotient_kind": quotient.kind,
        })
        input_times.clear()
        update_times.clear()

    try:
        adaptive.divide(product, outsider)
    except ValueError:
        pass
    else:
        raise AssertionError("division by the absent prime succeeded")
    return rows


def prime_cases(checkpoints, pure_limit, packed_limit, seed):
    primes, sympy_version = prime_sample(max(checkpoints) + 1, seed)
    outsider = primes.pop()
    pure_checkpoints = tuple(count for count in checkpoints if count <= pure_limit)
    packed_checkpoints = tuple(count for count in checkpoints if count <= packed_limit)
    results = {
        "id_bits": 129,
        "sample_seed": seed,
        "sympy_version": sympy_version,
        "sampling": "Seeded SymPy isprime rejection sampling of distinct 129-bit odd integers.",
    }
    if pure_checkpoints:
        pure_primes = primes[:max(pure_checkpoints)]
        results["pure"] = pure_prime_cases(
            pure_primes, outsider, pure_checkpoints, fast=False
        )
        results["pure_fast"] = pure_prime_cases(
            pure_primes, outsider, pure_checkpoints, fast=True
        )
    results["adaptive_without_mining"] = adaptive_prime_cases(
        primes, outsider, checkpoints, mining=False
    )
    if packed_checkpoints:
        results["adaptive_packed"] = adaptive_prime_cases(
            primes[:max(packed_checkpoints)], outsider, packed_checkpoints, mining=True
        )
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", nargs="+", type=int, default=DEFAULT_SEEDS)
    parser.add_argument("--config", nargs=2, type=int, action="append",
                        metavar=("DEPTH", "WIDTH"))
    parser.add_argument(
        "--prime-counts", nargs="+", type=int, default=(32, 128, 512, 2048, 8192)
    )
    parser.add_argument("--pure-prime-limit", type=int, default=128)
    parser.add_argument("--packed-prime-limit", type=int, default=32)
    parser.add_argument("--skip-deep", action="store_true")
    parser.add_argument("--skip-motif", action="store_true")
    parser.add_argument("--skip-primes", action="store_true")
    args = parser.parse_args()
    configs = args.config or DEFAULT_CONFIGS
    if any(depth < 1 or not 2 <= width <= 64 for depth, width in configs):
        parser.error("depth must be positive and width must be from 2 to 64")
    if min(args.prime_counts) < 2 or min(args.pure_prime_limit, args.packed_prime_limit) < 0:
        parser.error("prime checkpoints must be at least 2 and lane limits nonnegative")
    results = {
        "python": platform.python_version(),
        "machine": platform.machine(),
        "methods": {
            "timings": "Wall-clock microseconds with conversion, update, division and query timed separately.",
            "bytes": "Actual serialized DAG size including one record per reachable node.",
            "random_recursive_bases": "Positive recursive exponents with randomly signed terms and independently sampled root signs. Distinct positions are multiples of four. Shared templates use a separately seeded pool.",
            "oracle": "Independent raw expression evaluation using the observed powers-of-two residue cycle. Prime cases fold explicit fixture IDs modulo the query modulus.",
            "adaptive_eager": "Allow up to 1,000 signed normalization work units. Return the smaller of the eager and deferred results.",
            "prime_updates": "Sequential product of distinct sampled IDs. The unmined factor-tree lane runs every checkpoint. Stage medians include updates since the previous checkpoint.",
            "division": "Pure checked division validates divisibility by residue. Adaptive division uses a retained, already-packed exact factor. For the first seed in each deep configuration, canonical division recovers the pure quotient from only the two normalized input markings.",
        },
    }
    try:
        import numpy
    except ImportError:
        results["numpy_version"] = None
    else:
        results["numpy_version"] = numpy.__version__
    if not args.skip_deep:
        results["deep"] = deep_cases(configs, args.seed)
    if not args.skip_motif:
        results["motif_control"] = motif_control(args.seed[0])
        results["shared_prefix_control"] = shared_prefix_control(args.seed[0])
        results["independent_wide_control"] = independent_wide_control(args.seed[0])
    if not args.skip_primes:
        results["prime_products"] = prime_cases(
            sorted(set(args.prime_counts)), args.pure_prime_limit,
            args.packed_prime_limit, args.seed[-1]
        )
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
