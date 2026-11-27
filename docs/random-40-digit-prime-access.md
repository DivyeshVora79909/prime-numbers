# Access numbers built from random 40-digit primes

**All three trials passed.** Each trial gave Alice a different access list of 12,290 distinct primes drawn across the entire 40-digit range. The saved recursive exponent circuit represented the *exact* product of those IDs. Granting access, revoking access, reloading the changed number, and reversing each edit all worked.

**Randomness made little difference to the final size.** The three optimized labels occupied **6,055,316 to 6,067,065 bytes**. The earlier label, built from primes that shared their low 70 bits, occupied **6,066,846 bytes**. Two random labels were smaller and one was larger. All three took roughly **30 times** the space of their packed integer.

## One access check, from beginning to end

In trial A, Alice's prime ID was `4909305867477490416068453498166360224791`. Her second field was a single signed-power circuit for a **487,591-digit** access number. To ask whether Bob had access, the reader used his ID as the modulus. A zero remainder means his prime divides Alice's number.

| Operation on Alice's record | ID being checked | Exact remainder | Has access? |
| --- | --- | ---: | --- |
| Load original label and check Bob | `6297913873665372629309222046954044096129` | `0` | Yes |
| Divide by Bob's ID, encode the new label, and reload | Bob | `157127901495832741780909666598032616709` | No |
| Load original label and check Carol | `1540147111665682694523005535063805113443` | `909191532276090787054272718549055768699` | No |
| Multiply by Carol's ID, encode the new label, and reload | Carol | `0` | Yes |
| Check Alice's ID against her original label | Alice | `800912630152996136444118296006238606015` | No |

The grant and revocation each started from the original record. Granting Bob access again after revocation restored the original bytes. Revoking Carol after the grant did the same. Attempting exact division by Carol's ID in the original record raised an error. The reader performed multiplication, division and modulo on the merged signed representation. It had no factor list in the saved label.

## Storage across the three trials

Every row represents 12,290 members. “Canonical” is the project's recursive signed-power encoding before storage optimization. The lower bound rules out every smaller H2N1 encoding of that exact integer. The exact minimum lies between the bound and the file size we achieved.

| Trial (seed) | Decimal digits | Packed bytes | Canonical bytes | Optimized bytes | Lower bound | Maximum further saving |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Earlier structured IDs | 487,622 | 202,481 | 9,427,361 | 6,066,846 | 6,066,721 | 125 B |
| A (20261128) | 487,591 | 202,468 | 8,611,850 | 6,067,065 | 6,066,969 | **96 B** |
| B (20261129) | 487,620 | 202,480 | 8,334,101 | 6,063,949 | 6,063,871 | **78 B** |
| C (20261130) | 487,598 | 202,471 | 8,700,320 | 6,055,316 | 6,055,217 | **99 B** |

The optimizer saved **29.55%, 27.24% and 30.40%** from the respective canonical files. The earlier saving was **35.65%**. The random files already started smaller before optimization. Their canonical roots still needed **540,211**, **539,881** and **539,121** signed terms, close to the earlier **540,241**. Even a perfect H2N1 optimizer could save less than 100 additional bytes on each of the three original random labels.

The same optimizer ran on every input. Its two-state dynamic program prices a signed term using both the root reference and the cost of defining its exponent. It then rewrites the low prefix so commonly used exponent definitions get short references throughout the file. An odd low prefix can be rewritten this way for any pattern of its bits. The bound permits alternative definitions at all depths. The earlier report gives the derivation and the exact format assumptions. [Storage certificate and arithmetic kernel](prime-product-reachability-demo.md#why-the-size-certificate-covers-all-exponent-depths)

| Trial | Bytes after revoking one member | Bytes after granting one outsider | Reversing each edit |
| --- | ---: | ---: | --- |
| A | 6,063,967 | 6,062,307 | Original bytes restored twice |
| B | 6,065,213 | 6,070,616 | Original bytes restored twice |
| C | 6,053,821 | 6,065,435 | Original bytes restored twice |

Trial A's number grew when Carol was added, yet its optimized encoding shrank by 4,758 bytes. Multiplication changes the entire set of signed positions. Encoding size can move in either direction on an individual update.

## How the sample was generated

For each seed, Python drew odd integers uniformly from `10^39 + 1` through `10^40 − 1`. After a small-divisor filter, FLINT's `fmpz.is_prime()` performed a proven primality test. The sampler kept distinct accepted values in draw order. Every 40-digit prime was eligible. The fixed seeds make these benchmark draws reproducible. [FLINT primality test](https://flintlib.org/doc/fmpz/fmpz_is_prime.html)

The 12,292 entity IDs in each trial belonged to Alice, 12,290 members and one outsider. It also drew 32 extra prime IDs solely for negative checks. All **36,972** primes across the three fixtures were distinct. In each trial, all 12,292 entity IDs had distinct low-70-bit patterns. The earlier primes all had the same pattern.

| Trial | Odd candidates drawn for 12,324 distinct primes | Entity IDs equal to 1 / 3 modulo 4 |
| --- | ---: | ---: |
| A | 562,710 | 6,171 / 6,121 |
| B | 560,745 | 6,111 / 6,181 |
| C | 561,617 | 6,131 / 6,161 |

The UUID comparison depends on what is being counted. The prime-number-theorem integral across this interval estimates **9.85 × 10³⁷** prime choices, about **126.2 bits of choice**. A UUIDv4 has 122 random bits. UUIDv7 uses 48 bits for a millisecond timestamp and leaves 74 other variable bits, which may include counters or submillisecond precision. An individual 40-digit prime here takes **130–133 bits** to store as an ordinary integer, so it exceeds a 128-bit UUID field. [UUIDv4 and UUIDv7 layouts, RFC 9562](https://www.rfc-editor.org/rfc/rfc9562.html#section-5.4)

## How the complete numerical result was verified

Let `P` be the intended product and `N` the number decoded from the saved label. The verifier established **N = P** for all three labels without expanding either value into one large integer.

1. It checked divisibility by **every intended prime**. To save verifier time, it grouped at most 16 IDs per modulus and checked that each group gave remainder zero. The 769 checks per trial covered all 12,290 members. Since the primes are distinct, their full product `P` divides `N`.
2. Independently, it calculated the bit length of `P` by multiplying lower and upper bounds, rounding outward and retaining about 160 leading bits. It calculated the bit length and positivity of `N` from the decoded signed blocks using bounded carries. The two methods agreed at **1,619,741**, **1,619,836** and **1,619,764** bits respectively.
3. A positive multiple of `P` with the same bit length must be `P` itself. The next positive multiple, `2P`, already needs another bit. This also rules out extra copies of any member.
4. For Alice, the designated outsider and 32 additional absent IDs in each trial, the reader's nonzero remainders matched an independent modular fold over the intended members. These were **102 exact nonmember checks** in total.

The builder, reader, all-member verifier and stock-library verifier ran as separate processes. At 4, 8 and 16 IDs, each trial exercised the repository's public multiplication, storage optimizer, residue and exact division APIs. At full scale the experiment used the earlier report's convolution and bounded-block arithmetic. The repository's public decoder loaded each saved file, reproduced the builder's canonical SHA-256, and independently returned the correct member and outsider remainders.

## Query and update times

The three seeded jobs ran concurrently on the same machine. The ranges below describe those three runs. The reported direct query time for each target is the median of three calls. The prepared query time is the median of seven. Earlier benchmark timings were collected under different load, so the storage comparison is more informative than a small timing difference.

| Operation | Observed time |
| --- | ---: |
| Draw and prove 12,324 prime IDs, once per seed | 108.6–111.1 s |
| Merge 12,290 imported prime rows into the access number | 13.19–13.40 s |
| Optimize and encode the original number | 9.05–9.13 s |
| Resolve stored exponent definitions for direct queries, once | 2.93–3.12 s |
| Check one member or outsider directly on the serialized root | 0.951–1.037 s |
| Extract signed coordinates and prepare the query blocks, after loading | 0.834–0.902 s |
| Check one member or outsider with prepared blocks | **4.5–6.6 ms** |
| Divide by one member ID, then optimize and encode | 0.234–0.246 s + 8.85–8.99 s |
| Multiply by one outsider ID, then optimize and encode | 0.209–0.215 s + 8.34–8.39 s |
| Audit every member and all 34 nonmembers after preparing blocks | 17.34–17.38 s |
| Load the large file with the independent stock decoder | 102.0–102.8 s |

Every prepared query visited **6,328 blocks** of 256 signed positions. The prepared tuples used 430,312–430,328 extra bytes per label. For these inputs, randomizing the primes preserved exact access checks, while optimized writes still took seconds and a warm query still scanned the label.

## Reproduce the experiment

Run the command below from the repository root. It uses the unchanged arithmetic source blocks from the earlier report, checks their hashes, and extracts this report's sampler and trial driver. The command installs the pinned primality library in a temporary directory. The temporary fixtures and label files are removed at exit. The recorded run used repository commit `6bee6d2`, Python 3.12.8, python-flint 0.8.0 (FLINT 3.3.1), and GCC 13.3.0. Set `TRIAL_WORKERS=1` for serial runs on a smaller machine.

Each seed should end with `PASS`, `EXACT PRODUCT PASS` and `STOCK REFERENCE PASS`. The following SHA-256 values identify the optimized labels.

| Trial | SHA-256 |
| --- | --- |
| A | `784af4da1badce918ff0a4e6eed5e55d7600179001942ba245afd56d6fd8576a` |
| B | `23a3fc928e20c83761d14cccfd488c9836b1a3de84cefbe7e96c4e62bfc035ea` |
| C | `fb58e16829ffdf8d73aad11f32cf36f3312f933d36fdb6c88675e406584b7940` |

```sh
python3 - <<'PY'
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import os
from pathlib import Path
import re
import subprocess
import sys
from tempfile import TemporaryDirectory

repo = Path.cwd()
prior = (repo/'docs/prime-product-reachability-demo.md').read_text()
current = (repo/'docs/random-40-digit-prime-access.md').read_text()
seeds = (20261128, 20261129, 20261130)


def blocks(document, language):
    marker = chr(96) * 3
    return re.findall('^' + marker + language + r'\n(.*?)^' + marker + '$',
                      document, re.M | re.S)


with TemporaryDirectory() as tmp:
    work = Path(tmp)
    for lang, name, digest in [
        ('python', 'access.py', '5c3ef77282da8255442ede7fa5295e3475769325166316185c84132e596cfc79'),
        ('cpp', 'merge.cpp', 'c8cc752249c93ba049b82ba3ff855d9f77747f9eac1d17477abf84987c0036ec')
    ]:
        source, = blocks(prior, lang)
        assert sha256(source.encode()).hexdigest() == digest
        (work/name).write_text(source)
    sampler, trial = blocks(current, 'python')
    (work/'sample.py').write_text(sampler)
    (work/'trial.py').write_text(trial)
    deps = work/'deps'
    subprocess.run([sys.executable, '-m', 'pip', 'install', '--no-deps', '--only-binary=:all:',
                    '--target', str(deps), 'python-flint==0.8.0'], check=True)
    subprocess.run(['c++', '-O3', '-std=c++17', str(work/'merge.cpp'),
                    '-o', str(work/'merge')], check=True)
    env = dict(os.environ, PYTHONPATH=os.pathsep.join((str(deps), str(repo/'src'))),
               PYTHONDONTWRITEBYTECODE='1')

    def run(seed):
        commands = [[sys.executable, '-u', str(work/'sample.py'), str(seed)]]
        commands += [[sys.executable, '-u', str(work/'trial.py'), str(seed), phase]
                     for phase in ('build', 'read', 'audit', 'reference')]
        for command in commands:
            result = subprocess.run(command, env=env, check=True, text=True, stdout=subprocess.PIPE)
            print(f'[{seed}]', result.stdout, flush=True)

    with ThreadPoolExecutor(max_workers=int(os.environ.get('TRIAL_WORKERS', '3'))) as pool:
        list(pool.map(run, seeds))
PY
```

<details>
<summary>Sampler for unrestricted 40-digit primes</summary>

```python
from hashlib import sha256
from pathlib import Path
from random import Random
from time import perf_counter
import sys
from flint import fmpz

RECORDS = 12292
EXTRAS = 32
BASE = Path(__file__).resolve().parent


def sample(seed):
    rng = Random(seed)
    first_odd, stop = 10**39 + 1, 10**40
    sieve = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
    ids, seen = [], set()
    draws = 0
    start = perf_counter()
    while len(ids) < RECORDS + EXTRAS:
        p = rng.randrange(first_odd, stop, 2)
        draws += 1
        if p in seen or any(p % s == 0 for s in sieve):
            continue
        if fmpz(p).is_prime():
            ids.append(p)
            seen.add(p)
    raw = b''.join(p.to_bytes(17, 'big') for p in ids)
    (BASE / f'fixture-{seed}.bin').write_bytes(raw)
    entities = ids[:RECORDS]
    print('sample', seed, 'seconds', round(perf_counter()-start, 3), 'draws', draws,
          'fixture SHA256', sha256(raw).hexdigest(),
          'mod4 counts', [sum(p % 4 == r for p in entities) for r in (1, 3)],
          'low70 distinct', len({p % (1 << 70) for p in entities}), flush=True)


if __name__ == '__main__':
    sample(int(sys.argv[1]))
```

</details>

<details>
<summary>Trial driver and independent exact-product audit</summary>

```python
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR
from hashlib import sha256
from math import prod
from pathlib import Path
from time import perf_counter
import json
import sys

import access as a
from sample import BASE, EXTRAS, RECORDS


def fixture(seed):
    data = (BASE / f'fixture-{seed}.bin').read_bytes()
    assert len(data) == (RECORDS + EXTRAS) * 17
    ids = [int.from_bytes(data[i:i+17], 'big') for i in range(0, len(data), 17)]
    assert len(set(ids)) == len(ids) and all(10**39 <= p < 10**40 for p in ids)
    members = [ids[1], ids[2], *ids[4:RECORDS]]
    return ids, members


def product_widths(members):
    down = Context(prec=24, rounding=ROUND_FLOOR)
    up = Context(prec=24, rounding=ROUND_CEILING)
    decimal_lo = decimal_hi = Decimal(1)
    lo = hi = 1
    shift = 0
    for p in members:
        decimal_lo = down.multiply(decimal_lo, Decimal(p))
        decimal_hi = up.multiply(decimal_hi, Decimal(p))
        lo, hi = lo * p, hi * p
        trim = max(0, hi.bit_length() - 160)
        lo >>= trim
        hi = (hi + (1 << trim) - 1) >> trim
        shift += trim
    assert decimal_lo.adjusted() == decimal_hi.adjusted()
    assert lo.bit_length() == hi.bit_length()
    return 1 + decimal_lo.adjusted(), shift + lo.bit_length()


def build(seed):
    ids, members = fixture(seed)
    owner, member, outside = ids[0], ids[4100], ids[3]
    decimals, bits = product_widths(members)
    print('members', len(members), 'decimal digits', decimals, 'exact bits', bits,
          'packed bytes', (bits + 7)//8, flush=True)
    facts = dict(owner=owner, member=member, outside=outside, owner_r=a.fold(members, owner),
                 outside_r=a.fold(members, outside),
                 removed_r=a.fold([p for p in members if p != member], member))
    a.builtin_check(ids)
    c = a.Circuit()
    rows = a.timed('Import prime signed-power rows', lambda: [a.input_row(c, p) for p in members])
    row = a.timed('Build merged product', lambda: a.merge(rows))
    del rows
    canonical = a.timed('Encode canonical product', lambda: a.encode(row))
    facts['canonical_sha'] = sha256(canonical).hexdigest()
    wire, bound = a.timed('Optimize and encode product', lambda: a.optimize(row))
    print('positions', len(row), 'canonical terms', len(row)-row.count(b'0'),
          'canonical bytes', len(canonical), 'weighted bytes', len(wire), 'lower', bound)
    print('canonical SHA256', facts['canonical_sha'], 'weighted SHA256', sha256(wire).hexdigest())
    assert 0 <= len(wire) - bound <= 125  # width and header range of these fixtures
    (a.HERE/'access.h2n1').write_bytes(wire)
    (a.HERE/'checks.json').write_text(json.dumps(facts))


def audit(seed):
    ids, members = fixture(seed)
    blocks = tuple(a.signed_blocks(a.dense(a.split((a.HERE/'access.h2n1').read_bytes()))))
    # Recover only the bit length of the saved value. Carry and digit are bounded.
    radix = 1 << a.WIDTH
    carry = width = 0
    for i, x in enumerate(reversed(blocks)):
        carry, digit = divmod(x + carry, radix)
        if digit:
            width = i * a.WIDTH + digit.bit_length()
    assert carry >= 0
    if carry:
        width = len(blocks) * a.WIDTH + carry.bit_length()
    assert width > 0 and width == product_widths(members)[1]

    start = perf_counter()
    for i in range(0, len(members), 16):
        modulus = prod(members[i:i+16])
        assert a.mod_blocks(blocks, modulus) == 0
    controls = [ids[0], ids[3], *ids[RECORDS:]]
    for p in controls:
        assert a.mod_blocks(blocks, p) == a.fold(members, p) != 0
    print('EXACT PRODUCT PASS', len(members), 'members,', len(controls),
          'absent IDs,', width, 'bits, scan seconds', round(perf_counter()-start, 3))


def reference():
    facts = json.loads((a.HERE/'checks.json').read_text())
    c = a.Circuit()
    wire = (a.HERE/'access.h2n1').read_bytes()
    root = a.timed('Stock decode', lambda: a.decode_storage(c, wire))
    assert sha256(c.encode(root)).hexdigest() == facts['canonical_sha']
    assert a.residue(c, root, facts['member']) == 0
    assert a.residue(c, root, facts['outside']) == facts['outside_r']
    print('STOCK REFERENCE PASS')


if __name__ == '__main__':
    seed, mode = int(sys.argv[1]), sys.argv[2]
    a.HERE = BASE / str(seed)
    a.HERE.mkdir(exist_ok=True)
    if mode == 'build':
        import shutil
        shutil.copy2(BASE/'merge', a.HERE/'merge')
        build(seed)
    elif mode == 'read':
        a.read_and_update()
    elif mode == 'audit':
        audit(seed)
    elif mode == 'reference':
        reference()
    else:
        raise ValueError(mode)
```

</details>
