# Prime Numbers

Exact symbolic arithmetic for large integers represented as recursively signed powers of two. The project explores using prime products as membership labels and measures how different representations affect arithmetic, storage and queries.

- Canonical signed base-2 circuits support addition, subtraction, multiplication and exact division without materializing the represented integers.
- The grouped multiplication kernel reuses high parts of exponents and applies signed convolution to bounded offsets.
- Adaptive expressions can keep repeated blocks, powers, geometric sums and products factored. They support serialization, residue queries and cancellation of known factors.
- Seeded benchmarks cover independent deep inputs, shared exponent structure and products of sampled 129-bit prime IDs.

## Setup

Requires Python 3.11 or later.

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
python -m unittest discover -s tests -v
```

On Windows PowerShell, activate the environment with `.venv\Scripts\Activate.ps1`.

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

## Benchmarks

Install the benchmark extras to include the prime-product workload and the optional NumPy acceleration.

```sh
python -m pip install -e '.[benchmark]'
python -m symbolic_primes.deep_benchmark
python -m symbolic_primes.quotient_benchmark
```

See the [benchmark methods and results](docs/benchmarks.md) and [algorithm notes](docs/algorithms.txt).

## Scope

This is a research implementation. Canonical signed form minimizes the number of top-level signed terms. The adaptive search compares encoded sizes within a bounded set of rewrites. Retaining a product can speed updates and exact factor cancellation while increasing the cost of general residue queries. The decoder expects trusted output from the project encoder.
