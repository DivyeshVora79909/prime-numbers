# Benchmarks

Run from the repository root with Python 3.11 or later.

```sh
python -m pip install -e '.[benchmark]'
python -m symbolic_primes.benchmark
python -m symbolic_primes.deep_benchmark
python -m symbolic_primes.quotient_benchmark
```

Each benchmark checks its results and emits JSON. Products and quotients stay in symbolic representations throughout the runs. The archived runs used Python 3.12.14 on x86-64 with NumPy 2.3.5 and SymPy 1.14.0. Timings are wall-clock measurements. Storage figures count serialized bytes in a single reachable DAG. Fixture construction and verification run outside the operation timings. [Initial arithmetic](../benchmarks/results/benchmark-results.json) · [Deep arithmetic and prime products](../benchmarks/results/deep-benchmark-results.json) · [Exact quotient recovery](../benchmarks/results/quotient-benchmark-results.json)

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
