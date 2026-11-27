# One access number for a set of 12,290 prime IDs

**PASS for factor-free grants, revocations and membership checks.** An entity owns a prime ID and one access number. In this example Alice stores the product of the IDs allowed to access her. To check Bob's access, ask whether `Alice.access mod Bob.id == 0`. Each ID occurs once. Revoking Bob divides Alice's number by his ID. Granting a new ID multiplies it in after checking that it is absent.

| Alice's record | Value |
| --- | --- |
| Own prime ID | `6884203189328034036601239532666145996801` |
| Access number | One H2N1 recursive signed-power circuit, representing the product of 12,290 other distinct 40-digit IDs |
| Size of the represented number | 487,622 decimal digits, or 202,481 packed binary bytes |
| Stored circuit | 6,066,846 bytes, 540,262 signed terms in the root |

The earlier version of this report retained the original primes as leaves of a multiplication tree. Its divisibility results were correct. That choice left the merged representation untested. This run multiplied the signed rows together. The only objects in the saved number are shared definitions of the form `sum(sign × 2^exponent)`. Every exponent has a definition of the same form, down to zero.

## Results from the saved access number

The builder exited before the reader started. The reader had the stored circuit, three individual IDs and the expected remainders used by the assertions. The membership functions take only the representation and the queried ID. The original list of members was available to the builder for an independent modular check.

| Operation | Queried prime ID | Exact remainder | Access |
| --- | --- | ---: | --- |
| Read Alice's original number | `2635926831623976891135893021302819127297` | `0` | Yes |
| Revoke that member by exact division, save and reload | Same ID | `1949858375040160616871536411725593580240` | No |
| Read Alice's original number | `2106167842412050874497823083008505151489` | `247567763686552324532647455335898749986` | No |
| Grant that outsider by multiplication, save and reload | Same ID | `0` | Yes |
| Query Alice's own ID against her original number | `6884203189328034036601239532666145996801` | `488051596906232125867934281965766815158` | No |

The two edits started independently from Alice's original number. Reversing **each** edit restored all **6,066,846 original bytes** after reoptimization. Attempted exact division by the absent prime raised an error. The zero-remainder rule is exact because a prime divides a product of distinct primes precisely when it is one of them.

## Time and space

Measured on this machine with Python 3.12.8 and a C++17 multiplication kernel. Direct query times are medians of three calls. Prepared query times are medians of seven. Build and update times below are individual measurements.

| Operation | Measured time | What the operation reads or writes |
| --- | ---: | --- |
| Merge the 12,290 imported signed-power rows | 10.685 s | Fully combines intermediate products |
| Optimize and encode the original product | 8.119 s | Prices the whole label, verifies the rewrite and emits its bytes |
| Load exponent definitions once with the demonstration reader | 2.791 s | Resolves the small exponent coordinates, leaves the root serialized |
| Check member / outsider directly on the serialized root | 1.037 s / 1.020 s | Scans the signed references and accumulates modular powers |
| Prepare repeated queries after loading | 0.876 s + 0.030 s | Extracts the signed root coordinates and groups them into 6,328 blocks |
| Check member / outsider using the prepared blocks | 4.7 ms / 4.6 ms | Scans all blocks. The block tuple occupies 430,328 extra bytes on this Python runtime |
| Revoke one member | 0.256 s + 8.776 s | Exact arithmetic, then optimization and encoding |
| Grant one outsider | 0.170 s + 6.973 s | Merged arithmetic, then optimization and encoding |

Granting and revoking need seconds if each update must produce the optimized bytes immediately. A warm query is faster because the exponent coordinates have already been resolved. Its work still grows with the number of blocks. For these dense products at fixed ID width, the bit span grows roughly in proportion to the number of members. Eliminating the list of primes therefore gives no general guarantee of sublinear access checks.

| Representation of Alice's original number | Bytes |
| --- | ---: |
| Packed binary magnitude | 202,481 |
| Previous factor-retaining representation | 432,869 |
| Canonical factor-free recursive signed powers | 9,427,361 |
| Optimized factor-free recursive signed powers | **6,066,846** |
| Proven lower bound for every encoding in the same H2N1 format | **6,066,721** |

The optimized circuit saves 3,360,515 bytes over the canonical circuit. Even an exactly optimal H2N1 encoder would need about **30 times the packed size** for this particular number. The experiment certifies a gap of **at most 125 bytes**, or **0.0021%**, to the smallest valid H2N1 encoding. The exact minimum within that interval is unresolved. Packed magnitude counts only the integer payload, and the query cache is extra.

## How the operations work

- **Multiply.** A term at exponent `i` multiplied by a term at exponent `j` contributes their sign product at exponent `i+j`. The kernel combines equal positions using exact polynomial convolution, then propagates carries to yield one signed row. It uses direct convolution for a short operand and a number-theoretic transform for wider pairs. At the supported span of at most 2²³ convolution coefficients, each unreduced coefficient has magnitude below half the transform modulus 998244353, so its sign and exact integer value are recoverable. The output contains no multiply instruction or factor identity.
- **Divide.** Read the merged signed row from high positions to low in blocks of 256. If `r` is the previous remainder and `b` the signed value of the next block, `divmod(r × 2²⁵⁶ + b, p)` gives a quotient block and the next remainder. The final remainder must be zero. Normalizing the quotient blocks produces a new independent circuit. All tests here divide by an individual 40-digit ID.
- **Check.** A direct query computes the signed sum of powers modulo the query ID, advancing the power by the gap between exponents. A prepared query uses the same 256-position recurrence as division and keeps only the remainder. The prepared blocks depend on the circuit alone and can be reused for a previously unseen ID.
- **Verify.** The test generated 12,292 distinct certified primes with a fixed seed. Each has the Proth form `p = k × 2^70 + 1` with odd `k < 2^70` and witness `3^((p−1)/2) ≡ −1 (mod p)`. Thus the low 70 bits are shared across these IDs. The decimal width is certified by multiplying outward-rounded 24-digit bounds. The test never constructs the full product as a Python or C++ integer or decimal string. It does allocate arrays spanning the product's 1,619,845 signed-digit positions.

## Why the size certificate covers all exponent depths

Let `w(e)` be the minimum number of signed powers of two needed for exponent `e`, given by its nonadjacent form. If the root selects exponent `e`, H2N1 must store a reference to it, a row-length field for its definition, and at least `w(e)` child references in that definition. Distinct root positions require distinct rows. Each field or inner reference costs at least one byte even if other parts of the circuit reuse every child.

For a root with `m` distinct references, those references require at least `3m − 64 − 8192` bytes. Only 64 distinct row IDs can use one byte, and only 8,192 can use at most two. This gives an optimistic cost of `4 + w(e)` per selected position, less the fixed reference discounts. The root's minimum signed weight is 540,241, so the file header and the root-arity field together need at least 10 more bytes.

The repository's two-state carry dynamic program minimizes that weighted sum over signed root expansions in a linear pass. A position needs only the incoming carry, either zero or one. Considering positions through one above the supplied canonical row is sufficient. In a positive representation reaching beyond the integer's bit length, its two highest positions must have opposite signs on adjacent powers. Combining them eliminates the highest term and reduces the positive weighted cost. The minimum weighted sum in this run was 6,074,967. Subtracting 64 and 8,192 and adding 10 gives **6,066,721**. Allowing free sharing at deeper levels makes this a valid lower bound on every legal circuit.

To construct the upper bound, the encoder rewrites the low signed digits so that positions 0 through 21 occur early in the root. On this input, the low `+1` becomes `2²¹ − 2²⁰ − … − 2¹ − 2⁰`. This spends 21 extra root terms. It puts the definitions used throughout the rest of the circuit into cheap one-byte reference slots. The encoder measures the resulting complete file and verifies that normalizing it returns every original signed position.

## What runs in the repository and what runs in this report

The existing `multiply_storage → decode_storage → residue → divide_fast → optimize_storage` pipeline passed at 4, 8 and 16 IDs with a 192-state search budget. At each checkpoint, true division recovered the preceding canonical value and reoptimization reproduced the preceding stored bytes. At full scale, the repository's public `decode_storage` independently loaded the optimized output, `Circuit.encode` reproduced the canonical product's SHA-256, and `residue` reproduced the outsider remainder.

The full-size arithmetic test uses the report's exact convolution kernel, bounded-block division and ordinary-coordinate reader together with the repository's weighted carry DP. The large run exercises positive odd access numbers with ordinary exponent positions. Its results establish neither arbitrary-height performance nor an exactly optimal encoder for every integer. The application package was not modified.

## Reproduce

Run this command from the repository root with Python 3.11+ and a C++17 compiler called `c++`. The command extracts the two source blocks below into a temporary directory. It runs the builder and the access checks in separate processes, followed by the slower stock-decoder cross-check. The latter can take about 80–100 seconds and about 1 GiB of memory on this workload. The test prints `PASS` and `Full repository decoder and residue PASS`.

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 - <<'PY'
from pathlib import Path
from tempfile import TemporaryDirectory
import re
import subprocess
import sys

doc = Path("docs/prime-product-reachability-demo.md").read_text()
with TemporaryDirectory() as tmp:
    work = Path(tmp)
    for lang, name in [("python", "access.py"), ("cpp", "merge.cpp")]:
        block, = re.findall(r"^```" + lang + r"\n(.*?)^```$", doc, re.M | re.S)
        (work/name).write_text(block)
    subprocess.run(["c++", "-O3", "-std=c++17", str(work/"merge.cpp"), "-o", str(work/"merge")], check=True)
    for mode in ("build", "read", "reference"):
        subprocess.run([sys.executable, "-u", str(work/"access.py"), mode], check=True)
PY
```

<details>
<summary>Python driver, storage and access checks</summary>

```python
from array import array
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR
from hashlib import sha256
from pathlib import Path
from random import Random
from statistics import median
from time import perf_counter
import json
import subprocess
import sys

from symbolic_primes import Circuit, decode_storage, divide_fast, multiply_storage, optimize_storage, residue
from symbolic_primes._storage_bound import _naf_weights, _root_digits

HERE = Path(__file__).resolve().parent
WIDTH = 256
POS = bytes.maketrans(b'+-', b'10')
NEG = bytes.maketrans(b'+-', b'01')

def certified_primes(count):
    rng = Random(20261128)
    scale = 1 << 70
    first = (10**39 - 1 + scale - 1) // scale
    first += (1 - first) % 6
    last = (10**40 - 2) // scale
    result, seen = [], set()
    small_primes = (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
    while len(result) < count:
        k = rng.randrange(first, last + 1, 6)
        p = k * scale + 1
        if p in seen or any(p % d == 0 for d in small_primes):
            continue
        if pow(3, (p - 1) // 2, p) == p - 1:
            assert k % 2 == 1 and k < scale and 10**39 <= p < 10**40
            seen.add(p)
            result.append(p)
    return result


def timed(label, operation, repeats=1):
    elapsed = []
    for _ in range(repeats):
        start = perf_counter()
        result = operation()
        elapsed.append(perf_counter() - start)
    print(f'{label}: {median(elapsed):.4f} s', flush=True)
    return result


def merge(rows):
    payload = b'\n'.join(rows) + b'\n'
    return subprocess.run([HERE/'merge'], input=payload, capture_output=True, check=True).stdout


def input_row(circuit, p):
    terms = circuit.integer(p).terms
    row = bytearray(b'0') * (terms[-1][0]._small + 1)
    for exponent, sign in terms:
        row[exponent._small] = ord('+') if sign > 0 else ord('-')
    return bytes(row)


# H2N1 codec for this finite coordinate range. No row for a factor is retained.
def terms(n):
    i = 0
    while n:
        if n & 1:
            d = 2 - n % 4
            yield i, d
            n -= d
        n //= 2
        i += 1

def uvar(buf, x):
    while x >= 128:
        buf.append(128 | (x & 127))
        x >>= 7
    buf.append(x)

def encode(row):
    indices, body, top = {}, bytearray(), bytearray()
    def intern(e):
        if e not in indices:
            refs = [2*intern(i)+(d < 0) for i,d in terms(e)]
            indices[e] = len(indices)
            uvar(body, len(refs))
            for ref in refs: uvar(body, ref)
        return indices[e]
    uvar(top, len(row)-row.count(b'0'))
    for i, digit in enumerate(row):
        if digit != 48: uvar(top, 2*intern(i)+(digit == 45))
    header = bytearray(b'H2N1')
    uvar(header, len(indices)+1)
    return bytes(header+body+top)

def split(wire):
    assert wire[:4] == b'H2N1'
    view = memoryview(wire)
    offset = 4
    def read():
        nonlocal offset
        x = shift = 0
        while True:
            b = view[offset]
            offset += 1
            x |= (b & 127) << shift
            if b < 128: return x
            shift += 7
    nodes = array('I')
    count = read()
    for _ in range(count-1):
        value = 0
        for _ in range(read()):
            edge = read()
            e = nodes[edge//2]
            assert e <= 23
            value += (1 if edge % 2 == 0 else -1) << e
        assert 0 <= value <= 1 << 23
        nodes.append(value)
    n = read()
    return nodes, n, view[offset:]

def root_edges(parts):
    nodes, n, data = parts
    x = shift = 0
    for b in data:
        x |= (b & 127) << shift
        if b < 128:
            yield nodes[x//2], 1 if x%2 == 0 else -1
            x = shift = 0
            n -= 1
        else:
            shift += 7
    assert n == shift == 0

def dense(parts):
    row = bytearray()
    for e,s in root_edges(parts):
        assert e >= len(row)
        row.extend(b'0'*(e-len(row)))
        row.append(45 if s < 0 else 43)
    return bytes(row)


def optimize(row):
    # Canonical input, positive odd value, ordinary exponent positions.
    bits = bytearray(len(row) + 1)
    carry = 0
    for i, symbol in enumerate(row):
        carry += (symbol == 43) - (symbol == 45)
        bits[i] = carry % 2
        carry //= 2
    assert carry == 0
    costs = bytearray(4 + w for w in _naf_weights(len(row)))
    chosen = _root_digits(bits, costs)  # repository's weighted carry DP
    cost = sum(costs[e] for e, sign in chosen)
    best = bytearray(b'0') * (chosen[-1][0] + 1)
    for e, sign in chosen:
        best[e] = ord('+') if sign > 0 else ord('-')

    # Recode the EXACT low prefix with every position occupied. This moves
    # frequently used exponent definitions into cheap, early H2N1 row IDs.
    k = len(bits).bit_length()
    assert len(best) > k + 1 and bits[0] == 1
    prefix = sum(((s == 43) - (s == 45)) << i for i, s in enumerate(best[:k+1]))
    positive = (prefix + (1 << (k+1)) - 1) // 2
    for i in range(k + 1):
        best[i] = ord('+') if (positive >> i) & 1 else ord('-')
    assert merge([best, b'+']) == row
    wire = encode(best)
    root_weight = len(row) - row.count(b'0')
    field_bytes = lambda n: max(1, (n.bit_length() + 6) // 7)
    lower = 4 + field_bytes(root_weight) + field_bytes(root_weight+1) + cost - 64 - 8192
    assert lower <= len(wire)
    return wire, lower


def signed_blocks(row):
    for i in range((len(row)-1)//WIDTH*WIDTH, -1, -WIDTH):
        digits = row[i:i+WIDTH][::-1]
        yield int(digits.translate(POS), 2) - int(digits.translate(NEG), 2)


def mod_blocks(blocks, p):
    remainder = 0
    for block in blocks:
        remainder = ((remainder << WIDTH) + block) % p
    return remainder


def mod_wire(parts, p):
    power, previous, remainder = 1, 0, 0
    for exponent, sign in root_edges(parts):
        gap = exponent - previous
        power = ((power << gap) if gap <= WIDTH else power * pow(2, gap, p)) % p
        remainder = (remainder + sign * power) % p
        previous = exponent
    return remainder


def divide(row, p):
    assert p > 0
    quotient = bytearray(b'0') * len(row)
    remainder = 0
    for j, block in enumerate(signed_blocks(row)):
        digit, remainder = divmod((remainder << WIDTH) + block, p)
        assert abs(digit).bit_length() <= WIDTH
        offset = ((len(row)-1)//WIDTH - j) * WIDTH
        bits = format(abs(digit), 'b')[::-1].encode()
        quotient[offset:offset+len(bits)] = bits.replace(b'1', b'+' if digit >= 0 else b'-')
    if remainder:
        raise ValueError('not divisible')
    return merge([quotient, b'+'])  # normalize the quotient itself


def fold(members, p):
    r = 1
    for x in members:
        r = r * x % p
    return r


def builtin_check(ids):
    c = Circuit()
    wire = c.encode(c.one)
    for i, p in enumerate(ids[:16], 1):
        previous = wire
        result = multiply_storage(c, previous, c.integer(p), max_states=192)
        wire = result.data
        if i in (4, 8, 16):
            fresh = Circuit()
            root = decode_storage(fresh, wire)
            assert residue(fresh, root, ids[0]) == 0
            assert residue(fresh, root, ids[16]) == fold(ids[:i], ids[16]) != 0
            quotient = divide_fast(fresh, root, fresh.integer(p), max_terms=i*134)
            assert quotient is decode_storage(fresh, previous)
            assert optimize_storage(fresh, quotient, max_states=192).data == previous
            print('stock API', i, 'IDs,', len(wire), 'bytes, optimal', result.optimal)


def build():
    ids = timed('Generate certified prime IDs', lambda: certified_primes(12292))
    owner, member, outside = ids[0], ids[4100], ids[3]
    members = [ids[1], ids[2], *ids[4:]]
    assert len(set(members)) == len(members) == 12290
    assert owner not in members and member in members and outside not in members
    without = [v for v in members if v != member]
    facts = dict(owner=owner, member=member, outside=outside, owner_r=fold(members, owner),
                 outside_r=fold(members, outside), removed_r=fold(without, member))
    lo = hi = Decimal(1)
    down, up = Context(prec=24, rounding=ROUND_FLOOR), Context(prec=24, rounding=ROUND_CEILING)
    for p in members:
        lo, hi = down.multiply(lo, Decimal(p)), up.multiply(hi, Decimal(p))
    assert lo.adjusted() == hi.adjusted() == 487621
    builtin_check(ids)

    c = Circuit()
    rows = timed('Import prime signed-power rows', lambda: [input_row(c, x) for x in members])
    row = timed('Build merged product', lambda: merge(rows))
    del rows
    canonical = timed('Encode canonical product', lambda: encode(row))
    facts['canonical_sha'] = sha256(canonical).hexdigest()
    wire, bound = timed('Optimize and encode product', lambda: optimize(row))
    print('positions', len(row), 'canonical terms', len(row)-row.count(b'0'),
          'canonical bytes', len(canonical), 'weighted bytes', len(wire), 'global lower bound', bound)
    print('canonical SHA256', facts['canonical_sha'], 'weighted SHA256', sha256(wire).hexdigest())
    # Only the saved label and query/verification fixtures pass to the reader.
    (HERE/'access.h2n1').write_bytes(wire)
    (HERE/'checks.json').write_text(json.dumps(facts))


def read_and_update():
    facts = json.loads((HERE/'checks.json').read_text())
    member, outside, owner = (facts[k] for k in ('member', 'outside', 'owner'))
    wire = (HERE/'access.h2n1').read_bytes()
    parts = timed('Resolve stored exponents', lambda: split(wire))
    print('stored root terms', parts[1], 'wire nodes', len(parts[0])+1)
    for p, expected in ((member, 0), (outside, facts['outside_r']), (owner, facts['owner_r'])):
        r = timed('Direct wire query', lambda: mod_wire(parts, p), repeats=3)
        assert r == expected
        print('ID', p, 'remainder', r)
    row = timed('Extract signed coordinates', lambda: dense(parts))
    canonical = merge([row, b'+'])
    assert sha256(encode(canonical)).hexdigest() == facts['canonical_sha']
    blocks = timed('Prepare blocks', lambda: tuple(signed_blocks(row)))
    index_bytes = sys.getsizeof(blocks) + sum(sys.getsizeof(v) for v in {id(x):x for x in blocks}.values())
    print('prepared blocks', len(blocks), 'bytes', index_bytes)
    for p, expected in ((member, 0), (outside, facts['outside_r'])):
        assert timed('Prepared query', lambda: mod_blocks(blocks, p), repeats=7) == expected

    c = Circuit()
    for action, p in [('revoke', member), ('grant', outside)]:
        if action == 'revoke':
            changed = timed('Exact division', lambda: divide(row, p))
        else:
            assert mod_blocks(blocks, p) != 0
            changed = timed('Merged multiplication', lambda: merge([row, input_row(c, p)]))
        saved, lower = timed('Optimize and encode update', lambda: optimize(changed))
        updated = dense(split(saved))
        r = mod_blocks(signed_blocks(updated), p)
        print(action, len(saved), 'bytes, lower bound', lower, 'remainder', r)
        assert r == (facts['removed_r'] if action == 'revoke' else 0)
        restored = merge([updated, input_row(c, p)]) if action == 'revoke' else divide(updated, p)
        assert restored == canonical and optimize(restored)[0] == wire
    try:
        divide(row, outside)
    except ValueError:
        pass
    else:
        raise AssertionError('inexact revocation accepted')
    print('PASS')


if __name__ == '__main__':
    if sys.argv[1] == 'build':
        build()
    elif sys.argv[1] == 'read':
        read_and_update()
    elif sys.argv[1] == 'reference':
        facts = json.loads((HERE/'checks.json').read_text())
        c = Circuit()
        root = decode_storage(c, (HERE/'access.h2n1').read_bytes())
        assert sha256(c.encode(root)).hexdigest() == facts['canonical_sha']
        assert residue(c, root, facts['outside']) == facts['outside_r']
        print('Full repository decoder and residue PASS')
```

</details>

<details>
<summary>Exact signed-row multiplication kernel</summary>

```cpp
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
using Row = vector<int>;
constexpr int M = 998244353;  // 119 * 2^23 + 1; primitive root 3
int mul(int a, int b) { return int(int64_t(a) * b % M); }
int power(int a, int n) {
    int v = 1;
    while (n) { if (n & 1) v = mul(v, a); a = mul(a, a); n >>= 1; }
    return v;
}
void ntt(Row& a, bool inverse) {
    int n = a.size();
    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) swap(a[i], a[j]);
    }
    for (int w = 2; w <= n; w *= 2) {
        int step = power(3, (M-1)/w);
        if (inverse) step = power(step, M-2);
        for (int i = 0; i < n; i += w) {
            int t = 1;
            for (int j = 0; j < w/2; ++j) {
                int u = a[i+j], v = mul(a[i+j+w/2], t);
                a[i+j] = u+v < M ? u+v : u+v-M;
                a[i+j+w/2] = u-v >= 0 ? u-v : u-v+M;
                t = mul(t, step);
            }
        }
    }
    if (inverse) { int k = power(n, M-2); for (int& x : a) x = mul(x, k); }
}
Row normalize(Row v) {
    v.push_back(0);
    for (size_t i = 0; i+1 < v.size(); ++i) {
        int x = v[i], d = 0;
        if (x % 2) d = 2 - (((x + 2*v[i+1]) % 4 + 4) % 4);
        v[i] = d;
        v[i+1] += (x-d)/2;
        if (i+2 == v.size() && v.back()) v.push_back(0);
    }
    while (!v.empty() && !v.back()) v.pop_back();
    return v;
}
Row product(Row a, Row b) {
    if (a.empty() || b.empty()) return {};
    if (a.size() > b.size()) a.swap(b);
    int length = int(a.size()+b.size()-1);
    if (length > (1<<23)) throw runtime_error("coefficient span exceeds kernel limit");
    if (a.size() <= 140) {
        Row c(length);
        for (size_t i = 0; i < a.size(); ++i) if (a[i])
            for (size_t j = 0; j < b.size(); ++j) c[i+j] += a[i]*b[j];
        return normalize(move(c));
    }
    int n = 1;
    while (n < length) n *= 2;
    a.resize(n); b.resize(n);
    for (int& x : a) if (x < 0) x += M;
    for (int& x : b) if (x < 0) x += M;
    ntt(a, false); ntt(b, false);
    for (int i = 0; i < n; ++i) a[i] = mul(a[i], b[i]);
    ntt(a, true); a.resize(length);
    for (int& x : a) if (x > M/2) x -= M;
    return normalize(move(a));
}
int main() {
    vector<Row> rows;
    string line;
    while (getline(cin, line)) {
        Row r;
        for (char c : line) {
            if (c != '+' && c != '-' && c != '0') return 1;
            r.push_back(c == '+' ? 1 : c == '-' ? -1 : 0);
        }
        rows.push_back(move(r));
    }
    if (rows.empty()) return 1;
    while (rows.size() > 1) {
        vector<Row> next;
        for (size_t i = 0; i < rows.size(); i += 2)
            next.push_back(i+1 == rows.size() ? move(rows[i]) : product(move(rows[i]), move(rows[i+1])));
        rows.swap(next);
    }
    for (int d : rows[0]) cout.put(d > 0 ? '+' : d < 0 ? '-' : '0');
}
```

</details>
