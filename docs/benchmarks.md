# Benchmarks

Run from the repository root with Python 3.11 or later.

```sh
python -m pip install -e '.[benchmark]'
python -m symbolic_primes.benchmark
python -m symbolic_primes.deep_benchmark
python -m symbolic_primes.quotient_benchmark
```

Each benchmark checks its results and emits JSON. Products and quotients stay in symbolic representations throughout the runs. The archived runs used Python 3.12.14 on x86-64 with NumPy 2.3.5 and SymPy 1.14.0. Timings are wall-clock measurements. Storage figures count serialized bytes in a single reachable DAG. Fixture construction and verification run outside the operation timings. [Initial arithmetic](../benchmarks/results/benchmark-results.json) · [Deep arithmetic and prime products](../benchmarks/results/deep-benchmark-results.json) · [Exact quotient recovery](../benchmarks/results/quotient-benchmark-results.json)

## Updating a saved 30,000-prime product

The test starts with one H2N1 file representing the product of 30,000
distinct random 40-digit primes. Each call receives only the file and
one prime. The reference implementation is the archived Python driver
with its compiled C++ multiplier. The current Python API calls Rust in
the same process. [Checkpoint generator and input hashes](prime-access-storage-growth.md#measurement-and-verification)

| Operation, starting from saved bytes in memory | Previous Python/C++ driver | Python API with Rust | Speedup |
| --- | ---: | ---: | ---: |
| Check whether the last prime divides the product | 3.508 s | 0.102 s | 34.5× |
| Remove that prime and return an optimized H2N1 file | 11.703 s | 0.278 s | 42.0× |
| Multiply the saved quotient by the prime and return an optimized H2N1 file | 11.748 s | 0.254 s | 46.3× |

These are three-run wall-clock medians on an Intel Core i5-12450H,
x86-64 Linux, Python 3.12.8 and a release build with Rust 1.83.0. Each
pair ran in separate fresh processes pinned to the same logical CPU,
alternating which implementation ran first.
Both columns include reading the exponent definitions from bytes on every
call. The update columns also include arithmetic, optimization and
serialization. Disk reads, prime generation and the final hash comparisons
are outside the timers. The historical driver used the Python package at
commit `da3ca37`. Its timings include the original carry DP, subprocess
calls and internal normalization checks where applicable. These timings measure
the complete workflows. [Historical driver and C++ source](prime-product-reachability-demo.md)

The member remainder was zero. The outsider remainder was
`833120395990529729763735569905656999601` (0.102 s median). Removing
the member gave a 15,681,645-byte quotient and a nonzero member remainder.
Restoring it reproduced all **15,686,193 bytes** of the original file.
As an independent construction check, the test converted the 30,000
individual IDs to signed rows and combined pairs with the Rust multiplier.
All 29,999 multiplications ran to completion. The final row matched all
3,953,815 positions of the archived canonical row, and encoding it
reproduced the same file byte for byte. Each intermediate held the
combined signed digits of its product.
The lower bound was **15,686,102 bytes**, leaving at most 91 bytes of
possible improvement in H2N1. For comparison, packed binary needs
494,227 bytes for this value. [Storage measurements and certificate](prime-access-storage-growth.md)

To repeat the checks, build the native extension and use `wire` from the
checkpoint generator linked above.

```python
from hashlib import sha256
from statistics import median
from timeit import repeat
from symbolic_primes import Circuit, divide_storage, multiply_storage, residue

c = Circuit()
member = 3589358312608496142804550037657486194871
outsider = 1540147111665682694523005535063805113443
factor = c.integer(member)
query = lambda: residue(c, wire, member)
remove = lambda: divide_storage(c, wire, member, max_terms=4_000_000)
quotient = remove().data
restore = lambda: multiply_storage(c, quotient, factor)

assert query() == 0
assert residue(c, wire, outsider) == 833120395990529729763735569905656999601
assert residue(c, quotient, member) != 0
assert sha256(quotient).hexdigest() == '1b0106e9a999efdc8be37cb5c00e170c4080e0c97a6a22d0978423a59bae00b5'
result = restore()
assert result.data == wire and result.max_extra_bytes == 91
for name, operation in [('member', query), ('remove', remove), ('restore', restore)]:
    print(name, median(repeat(operation, number=1, repeat=3)))
```

For the previous column, use the archived `access.py` and compiled
`merge` from the historical driver with that checkout. With the signed
prime row prepared before timing, the corresponding calls were
`access.mod_wire(access.split(wire), member)`,
`access.optimize(access.divide(access.dense(access.split(wire)), member))`,
and `access.optimize(access.merge([access.dense(access.split(quotient)), prime_row]))`.

### Where the time went

The following measurements isolate operations that have different costs
from a full saved-file update.

| Input and operation | Before | With Rust | What was timed |
| --- | ---: | ---: | --- |
| 12,290-prime signed row × next 40-digit prime | 73.1 ms | 10.6 ms | In-process C++ and Rust multiplication, including returning the signed row |
| Two 810,000-position signed rows | 737.7 ms | 95.6 ms | Same kernel comparison, with byte-identical output |
| Weighted carry DP on 1,619,742 positions | 1.547 s | 0.244 s | Recovering 539,627 signed terms, including Python tuples |
| Default storage search on a 2,048-bit ordinary integer | 120 ms | 54 ms | Whole optimizer, identical bytes, bound and 445 search states |
| Default storage search on a 4,096-bit ordinary integer | 234 ms | 130 ms | Whole optimizer, identical bytes, bound and 555 search states |

The multiplication rows use seven-call medians. The larger pair takes
the first and last 810,000 positions of the 12,290-prime row. The DP
uses `Random(103).randbytes(span)` with each byte reduced modulo two
and position cost `4 + NAF_weight(position)`. The storage runs use
`Random(103).getrandbits(width) | (1 << (width - 1)) | 1` and three-call
medians. Deep sparse inputs keep the recursive algorithms. A modular
query on the saved dense product still reads the encoded graph, and an
unlimited minimum-size proof can still require combinatorial search.

For dense addition on fresh 100,000-position circuits, the five-run median
fell from 342 ms to 259 ms and subtraction fell from 338 ms to 268 ms.
These timings include constructing the Python result graph. The native
addition loop alone took 1.46 ms on the same input width.

## Recursive multiplication

The deep suite uses 36 seeded signed input pairs at depths 3, 5, and 8. It varies width and draws exponents from independent branches or a shared template pool. It checks input and product residues using a separate oracle that finds powers-of-two residue cycles. Both signed multiplication kernels produce byte-identical canonical results.

| Fixture | Original signed | Optimized signed | Deferred factor tree | Signed output |
| --- | ---: | ---: | ---: | ---: |
| Depth 3, width 8, independent branches | 109.21 ms | 7.39 ms | 10.40 µs | 4,265 bytes |
| Depth 8, width 2, independent branches | 20.79 ms | 4.62 ms | 9.89 µs | 2,767 bytes |
| 128 × 128 terms, shared high parts | 4.75 s | 18.46 ms | | 46,599 bytes |
| 128 × 128 terms, independent branches | 2.44 s | 245.98 ms | | 209,434 bytes |

The first two rows report medians over three seeds. Each wide row reports one seeded trial. The shared-high-part fixture groups 16,384 local coefficient pairs under 16 high-part pairs and produces 2,444 signed terms. The independent wide fixture retains all 16,384 signed terms. The factor-tree timings measure constructing the product expression while preserving the operands.

In a separate control, motif packing reduced the adaptive encoding of an injected repeated block from 362 to 138 bytes. Placing that block inside an exponent reduced the encoding from 365 to 141 bytes. In the initial tower suite, a difference of squares yielded two top-level terms and 26 bytes, while multiplying ten factors with separated supports yielded 1,024 top-level terms. [Initial arithmetic](../benchmarks/results/benchmark-results.json) · [Deep arithmetic](../benchmarks/results/deep-benchmark-results.json)

## Exact quotient recovery

Each trial copies only the signed numerator and divisor into a fresh arena. It computes a canonical signed quotient, then compares the complete encoding with the independently constructed expected result. The table gives medians of three trials. [Exact quotient recovery](../benchmarks/results/quotient-benchmark-results.json)

| Numerator terms | Divisor terms | Divisor | Quotient terms | Median division |
| ---: | ---: | --- | ---: | ---: |
| 48 | 6 | Negative odd | 8 | 15.74 ms |
| 128 | 8 | Negative even | 16 | 44.41 ms |
| 320 | 10 | Negative odd | 32 | 99.67 ms |

## Prime-product updates

The workload adds distinct 129-bit IDs sampled with SymPy's `isprime` predicate. It checks member and outsider residues, known-factor removal, and reinsertion. Reference checks fold the sampled factors modulo the query. Update medians summarize the additions since the previous checkpoint and exclude input conversion. Division and query columns are checkpoint timings. [Prime-product results](../benchmarks/results/deep-benchmark-results.json)

| IDs | Representation | Serialized DAG | Packed binary size | Median add | Exact division | Outsider residue |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 128 | Original signed | 55,022 B | 2,058 B | 59.98 ms | 8.98 ms | 9.76 ms |
| 128 | Optimized signed | 55,022 B | 2,058 B | 45.32 ms | 9.26 ms | 9.02 ms |
| 128 | Factor tree | 10,063 B | 2,058 B | 18.34 µs | 27.31 µs | 1.78 ms |
| 8,192 | Factor tree | 610,652 B | ≤132,096 B | 33.64 µs | 197.68 µs | 122.38 ms |

The packed-binary size is inferred from the leading digits at 128 IDs and bounded by the total input bit lengths at 8,192 IDs. The factor tree divides by removing a retained factor. The signed lanes check divisibility and retain a quotient expression. At 8,192 IDs, the median factor conversion took 260.51 µs and the explicit factor-list reference query took 1.53 ms. The tree makes updates cheaper by retaining factors, which it then evaluates during general residue queries.

## Bounded recursive storage

A separate local run used the default 2,048-state budget and 1% target.
With Python 3.12.8 and one `Random(20261128)` generator, it sampled eight
ordinary odd integers at each width in ascending order using
`rng.getrandbits(w) | (1 << (w - 1)) | 1`. Each value started in a fresh
circuit. Timings cover optimization, with conversion and round-trip
verification outside the timer. Savings compare against the canonical
recursive encoding. The gap is an upper bound on the remaining possible
savings in the same H2N1 format.

| Input bits | Median optimization | Median byte savings | Largest certified gap in bytes | Inputs within target |
| ---: | ---: | ---: | ---: | ---: |
| 128 | 2.9 ms | 2.20% | 3 | 6/8 |
| 256 | 5.0 ms | 1.84% | 3 | 8/8 |
| 512 | 14.2 ms | 6.37% | 6 | 8/8 |
| 1,024 | 61.4 ms | 10.12% | 7 | 8/8 |
| 2,048 | 101.3 ms | 9.57% | 3 | 8/8 |
| 4,096 | 183.5 ms | 12.91% | 7 | 8/8 |
| 8,192 | 337.7 ms | 12.53% | 10 | 8/8 |

All 56 encodings shrank and decoded exactly. The solver proved 21
minima. The two 128-bit inputs that reached the state limit had a
largest certified relative gap of about 1.36%. In the separate mixed
tower fixture with six repeated blocks and three irregular terms, the
memoized solver proved the 148-byte minimum in 536 states.
