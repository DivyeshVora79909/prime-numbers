# Prime Numbers

Exact symbolic arithmetic for large integers represented as recursively signed powers of two. The project explores using prime products as membership labels and measures how different representations affect arithmetic, storage and queries.

- Canonical signed base-2 circuits support addition, subtraction, multiplication and exact division without materializing the represented integers.
- The grouped multiplication kernel reuses high parts of exponents and applies signed convolution to bounded offsets.
- Adaptive expressions can keep repeated blocks, powers, geometric sums and products factored. They support serialization, residue queries and cancellation of known factors.
- Seeded benchmarks cover independent deep inputs, shared exponent structure and products of sampled 129-bit prime IDs.

## Setup

Requires Python 3.11 or later. Install a stable Rust toolchain to build the
optional native extension from source. The standard install builds it in
release mode when a compiler is available. Python remains usable without it.

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
python -m unittest discover -s tests -v
```

On Windows PowerShell, activate the environment with `.venv\Scripts\Activate.ps1`.

To check that a build includes Rust, run `python -c "from symbolic_primes import _native"`.
Run the Rust unit tests with `cargo test --locked --manifest-path native/Cargo.toml`.

## Example

```python
from symbolic_primes import Adaptive, Circuit, divide_fast, multiply_fast

circuit = Circuit()
exponent = circuit.power(circuit.integer(100))
power = circuit.power(exponent)
left = circuit.subtract(power, circuit.one)
right = circuit.add(power, circuit.one)

product = multiply_fast(circuit, left, right)
assert divide_fast(circuit, product, right) is left

adaptive = Adaptive(circuit)
packed = adaptive.pack(product)
restored = adaptive.decode(adaptive.encode(packed))
assert adaptive.residue(restored, 97) == adaptive.residue(packed, 97)
```

## Storage optimization

`optimize_storage(circuit, value)` starts with the canonical recursive
encoding and keeps whole-graph byte savings. The deterministic default
allows 2,048 search states and aims to stop within 1% of a proved lower
bound. If the budget expires first, the result still decodes exactly.
`data` holds the bytes, `max_extra_bytes` bounds how much smaller any
valid encoding could be, and `optimal` reports a proved minimum.
`decode_storage(circuit, data)` restores the value for arithmetic. Use
`target_gap_percent=0` to pursue an exact proof within a budget, or
`max_states=None` for an unlimited exact search. [Algorithm and proof scope](docs/algorithms.txt#L84)

`multiply_storage` and `divide_storage` accept saved H2N1 bytes or canonical
operands from the given circuit. Division also accepts an integer divisor.
Both feed the arithmetic result through the same optimizer used for a
directly constructed value. `residue` can query the saved file in memory.
Each operation needs only its operands.

```python
from symbolic_primes import Circuit, divide_storage, multiply_storage, optimize_storage, residue

c = Circuit()
prime_id = 3589358312608496142804550037657486194871
original = optimize_storage(c, c.integer(101))
added = multiply_storage(c, original.data, c.integer(prime_id))
assert residue(c, added.data, prime_id) == 0
removed = divide_storage(c, added.data, prime_id)
assert removed.data == original.data
assert residue(c, removed.data, prime_id) == 101
```

`divide_storage` shares the `divide_fast` quotient budget of 256 signed
digits by default. Pass a larger `max_terms` for a large quotient.
`QuotientWorkLimit` means the digit budget ran out. An established inexact
division raises `ValueError`.

## Benchmarks

Install the benchmark extras for the prime-product workload and the optional NumPy convolution fallback.

```sh
python -m pip install -e '.[benchmark]'
python -m symbolic_primes.deep_benchmark
python -m symbolic_primes.quotient_benchmark
```

See the [30,000-prime Rust benchmark](docs/benchmarks.md#updating-a-saved-30000-prime-product) and [algorithm notes](docs/algorithms.txt).

## Scope

This is a research implementation. Canonical signed form minimizes the number of top-level signed terms. The adaptive search compares encoded sizes within a bounded set of rewrites. Retaining a product can speed updates and exact factor cancellation while increasing the cost of general residue queries. The decoder expects trusted output from the project encoder.
