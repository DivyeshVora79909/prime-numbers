"""Reproducible arithmetic benchmarks with no expanded products or quotients."""

import argparse
import json
import platform
from random import Random
from statistics import median
from time import perf_counter_ns

from .symbolic import Circuit, exact_divide, residue


def measure(operation):
    """Return the operation's result and elapsed microseconds."""
    start = perf_counter_ns()
    result = operation()
    return result, (perf_counter_ns() - start) / 1_000


def prime_product_mod(factors, modulus):
    result = 1 % modulus
    for factor in factors:
        result = (result * factor) % modulus
    return result


def tower_exponent(circuit):
    return circuit.power(circuit.power(circuit.integer(64)))


def recursive_cases(repeats):
    multiply_times = []
    shift_times = []
    geometric_times = []
    residue_times = []
    for _ in range(repeats):
        circuit = Circuit()
        one = circuit.integer(1)
        exponent = tower_exponent(circuit)
        large_power = circuit.power(exponent)
        less = circuit.subtract(large_power, one)
        more = circuit.add(large_power, one)

        product, elapsed = measure(lambda: circuit.multiply(less, more))
        multiply_times.append(elapsed)
        twice_exponent = circuit.add(exponent, exponent)
        assert product == circuit.subtract(circuit.power(twice_exponent), one)
        assert twice_exponent == circuit.power(
            circuit.add(circuit.power(circuit.integer(64)), one)
        )

        shifted = circuit.multiply(product, large_power)
        shifted_quotient, elapsed = measure(
            lambda: circuit.divide_power(shifted, exponent)
        )
        shift_times.append(elapsed)
        assert shifted_quotient == product

        # Fermat gives an independent reference from a bounded exponent.
        modulus = 97
        exponent_mod = pow(2, 1 << 64, modulus - 1)
        expected = (pow(2, 2 * exponent_mod, modulus) - 1) % modulus
        assert residue(circuit, product, modulus) == expected

        quotient, elapsed = measure(lambda: exact_divide(circuit, product, 3))
        geometric_times.append(elapsed)
        observed, elapsed = measure(lambda: residue(circuit, quotient, modulus))
        residue_times.append(elapsed)
        assert observed == expected * pow(3, -1, modulus) % modulus
        # The quotient is a sum of E terms all congruent to 1 modulo 3.
        assert residue(circuit, quotient, 3) == 1
        assert type(quotient).__name__ == "Geometric"

    return {
        "E": "2^(2^64)",
        "repeats_with_fresh_circuit": repeats,
        "difference_of_squares": {
            "input_pair_terms": len(less.terms) * len(more.terms),
            "output": "2^(2E)-1",
            "expanded_bit_length": "2E",
            "output_stats": circuit.stats(product),
            "median_us": median(multiply_times),
        },
        "division_by_power_of_two": {
            "operation": "((2^(2E)-1)*2^E)/2^E",
            "input_stats": circuit.stats(shifted),
            "output_stats": circuit.stats(shifted_quotient),
            "median_us": median(shift_times),
        },
        "mersenne_division_by_three": {
            "operation": "(2^(2E)-1)/3",
            "minimum_top_level_signed_terms": "E",
            "result_type": type(quotient).__name__,
            "geometric_count_stats": circuit.stats(quotient.count),
            "median_divide_us": median(geometric_times),
            "median_mod_97_us": median(residue_times),
        },
    }


def separated_terms(levels):
    circuit = Circuit()
    one = circuit.integer(1)
    product = one
    modulus = 97
    expected = 1
    exponent = tower_exponent(circuit)
    exponent_mod = pow(2, 1 << 64, modulus - 1)
    rows = []
    for level in range(1, levels + 1):
        multiplier = circuit.add(circuit.power(exponent), one)
        product, elapsed = measure(lambda: circuit.multiply(product, multiplier))
        expected = expected * (1 + pow(2, exponent_mod, modulus)) % modulus
        assert len(product.terms) == 1 << level
        assert residue(circuit, product, modulus) == expected
        if level % 2 == 0 or level == levels:
            rows.append(
                {
                    "factors": level,
                    "output_stats": circuit.stats(product),
                    "multiply_us": elapsed,
                }
            )
        exponent = circuit.add(exponent, exponent)
        exponent_mod = (2 * exponent_mod) % (modulus - 1)
    return {
        "factor_at_index_i": "1 + 2^(E*2^i), with i starting at 0",
        "E": "2^(2^64)",
        "checkpoints": rows,
    }


def prime_sample(count, seed):
    try:
        from sympy import __version__, isprime
    except ImportError as error:
        raise SystemExit(
            "Install SymPy for the prime workload, or pass --skip-primes."
        ) from error

    rng = Random(seed)
    sampled = set()
    primes = []
    while len(primes) < count:
        candidate = (1 << 128) | rng.getrandbits(128) | 1
        if candidate not in sampled and isprime(candidate):
            sampled.add(candidate)
            primes.append(candidate)
    return primes, __version__


def prime_cases(checkpoints, seed, query_repeats):
    primes, sympy_version = prime_sample(max(checkpoints) + 1, seed)
    outsider = primes.pop()
    circuit = Circuit()
    product = circuit.integer(1)
    interval_times = []
    rows = []
    for count, prime in enumerate(primes, start=1):
        prime_symbols = circuit.integer(prime)
        pair_terms = len(product.terms) * len(prime_symbols.terms)
        product, elapsed = measure(lambda: circuit.multiply(product, prime_symbols))
        interval_times.append(elapsed)
        if count not in checkpoints:
            continue

        factors = primes[:count]
        query_times = []
        factor_fold_times = []
        for _ in range(query_repeats):
            expected, elapsed_fold = measure(
                lambda: prime_product_mod(factors, outsider)
            )
            observed, elapsed_query = measure(
                lambda: residue(circuit, product, outsider)
            )
            assert observed == expected and observed != 0
            query_times.append(elapsed_query)
            factor_fold_times.append(elapsed_fold)
        assert residue(circuit, product, prime) == 0
        highest_position = product.terms[-1][0]._small
        assert highest_position is not None
        # The sign of the first lower digit determines which side of 2**h
        # the product lies on. Its lower nonadjacent digits cannot reverse it.
        below_highest_power = len(product.terms) > 1 and product.terms[-2][1] < 0
        bit_length = highest_position if below_highest_power else highest_position + 1
        rows.append(
            {
                "prime_factors": count,
                "expanded_bit_length": bit_length,
                "packed_binary_bytes": (bit_length + 7) // 8,
                "expanded_bytes_upper_bound": (129 * count + 7) // 8,
                "output_stats": circuit.stats(product),
                "last_update_pair_terms": pair_terms,
                "last_multiply_us": elapsed,
                "interval_median_multiply_us": median(interval_times),
                "median_nonmember_residue_us": median(query_times),
                "median_explicit_factor_list_residue_us": median(factor_fold_times),
            }
        )
        interval_times.clear()

    removed = primes[len(primes) // 2]
    removed_symbols = circuit.integer(removed)
    remaining = [prime for prime in primes if prime != removed]
    division_times = []
    removed_query_times = []
    restore_times = []
    for _ in range(query_repeats):
        quotient, division_time = measure(
            lambda: exact_divide(circuit, product, removed)
        )
        removed_residue, query_time = measure(
            lambda: residue(circuit, quotient, removed)
        )
        restored, restore_time = measure(
            lambda: circuit.multiply(quotient, removed_symbols)
        )
        assert restored is product
        division_times.append(division_time)
        removed_query_times.append(query_time)
        restore_times.append(restore_time)
    assert removed_residue == prime_product_mod(remaining, removed)
    assert removed_residue != 0
    assert residue(circuit, quotient, outsider) == prime_product_mod(
        remaining, outsider
    )
    if remaining:
        assert residue(circuit, quotient, remaining[0]) == 0

    return {
        "prime_bit_length": 129,
        "seed": seed,
        "sympy_version": sympy_version,
        "prime_sampling": "Seeded rejection sampling with SymPy isprime for 129-bit odd integers.",
        "query_repeats": query_repeats,
        "checkpoints": rows,
        "removing_a_known_factor_from_the_signed_product": {
            "result_type": type(quotient).__name__,
            "median_divide_us": median(division_times),
            "median_query_removed_prime_us": median(removed_query_times),
            "median_restore_by_multiplication_us": median(restore_times),
            "queried_prime_is_no_longer_a_divisor": True,
            "query_requires_no_modular_inverse_of_divisor": True,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--prime-counts", nargs="+", type=int, default=[1, 2, 4, 8, 16, 32, 64, 128]
    )
    parser.add_argument("--seed", type=int, default=20261128)
    parser.add_argument("--repeats", type=int, default=21)
    parser.add_argument("--growth-levels", type=int, default=10)
    parser.add_argument("--skip-primes", action="store_true")
    args = parser.parse_args()
    if args.repeats < 1 or args.growth_levels < 1 or min(args.prime_counts) < 1:
        parser.error("counts and repeats must be positive")

    results = {
        "python": platform.python_version(),
        "machine": platform.machine(),
        "methods": {
            "timings": "Wall-clock microseconds. Excludes input construction, size measurement, and reference checks.",
            "multiply": "One sequential workload. Stage medians include updates since the previous checkpoint.",
            "residues": "Each residue call uses a fresh per-call cache.",
            "bytes": "Actual encoded signed DAG length. Packed size is inferred from the leading signed digits. The geometric entry measures the signed count payload.",
        },
        "recursive_cases": recursive_cases(args.repeats),
        "unavoidable_output_growth": separated_terms(args.growth_levels),
    }
    if not args.skip_primes:
        results["seeded_prime_product"] = prime_cases(
            set(args.prime_counts), args.seed, min(args.repeats, 9)
        )
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
