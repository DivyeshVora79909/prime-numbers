"""Benchmark exact quotient recovery from unrelated recursive markings."""

from argparse import ArgumentParser
from json import dumps
from random import Random
from statistics import median
from time import perf_counter_ns

from .quotient_kernel import divide_fast
from .similarity_kernel import multiply_fast
from .symbolic import Circuit


SEED = 10912


def _fixture(circuit, random, anchors, size, odd, sign):
    terms = {circuit.zero if odd else circuit.integer(2): 1}
    while len(terms) < size:
        first, second = sorted(random.sample(anchors, 2))
        anchor = circuit.add(second, first) if random.randrange(2) else circuit.subtract(second, first)
        position = circuit.add(anchor, circuit.integer(4 * random.randrange(-18, 19)))
        terms[position] = random.choice((-1, 1))
    ordered = sorted(terms.items())
    ordered[-1] = ordered[-1][0], sign
    return circuit._intern(tuple(ordered))


def _copy_to(circuit, *markings):
    copies = {}

    def copy(marking):
        if marking not in copies:
            copies[marking] = circuit._intern(tuple((copy(p), d) for p, d in marking.terms))
        return copies[marking]

    return tuple(copy(marking) for marking in markings)


def run(repeats=3):
    random = Random(SEED)
    source = Circuit()
    tower = source.integer(23)
    for _ in range(5):
        tower = source.power(tower)
    anchors = [source.power(source.add(tower, source.integer(3 * i))) for i in range(19)]

    cases = []
    for quotient_terms, divisor_terms, odd in ((8, 6, True), (16, 8, False), (32, 10, True)):
        expected = _fixture(source, random, anchors, quotient_terms, True, 1)
        divisor = _fixture(source, random, anchors, divisor_terms, odd, -1)
        numerator = multiply_fast(source, expected, divisor)
        expected_wire = source.encode(expected)
        samples = []
        for _ in range(repeats):
            circuit = Circuit()
            numerator_copy, divisor_copy = _copy_to(circuit, numerator, divisor)
            start = perf_counter_ns()
            quotient = divide_fast(circuit, numerator_copy, divisor_copy, max_terms=quotient_terms)
            samples.append((perf_counter_ns() - start) / 1_000_000)
            if circuit.encode(quotient) != expected_wire:
                raise AssertionError("division returned a different canonical quotient")
        cases.append({
            "numerator": source.stats(numerator),
            "divisor": source.stats(divisor),
            "divisor_parity": "odd" if odd else "even",
            "divisor_sign": -1,
            "quotient": source.stats(expected),
            "division_median_ms": median(samples),
            "division_samples_ms": samples,
        })
    return {"seed": SEED, "repeats": repeats, "cases": cases}


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    print(dumps(run(args.repeats), indent=2))


if __name__ == "__main__":
    main()
