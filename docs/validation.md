# File-backed validation

The release CLI rebuilt the product of 30,000 distinct random 40-digit primes twice. The first run used the saved member stream. The second regenerated the same stream from its seed. Importing, proving and granting the entire batch took 74.70 s and 132.44 s respectively, excluding input generation. Both labels occupied **494,227 bytes** and matched the historical ordinary-binary SHA-256 `b4a6c3bb41f08ef79c63ab6cbc0dcaa65db2405c6ace6819859366ec80021c5b`. The last ID had remainder zero, and the outsider returned `833120395990529729763735569905656999601`. Removing the member reproduced the independent quotient hash. Granting it back restored the complete original file. Repeated grants were idempotent.

## Two GiB on local storage

The second fixture was an actual 2 GiB file containing a deterministically generated dense positive cofactor. It measures arithmetic and I/O at scale. The cofactor's other factors have unknown primality and multiplicity. The harness checked that the proposed ID was absent, granted it, checked membership, revoked it and verified a byte-identical restoration against streamed SHA-256 and independent bounded-block Python arithmetic. The release worker had a **96 MiB process address-space limit** throughout.

| Public operation | Warm input, runs 1 / 2 | After requesting input eviction, runs 1 / 2 |
| --- | ---: | ---: |
| Check membership | 1.96 / 3.74 s | 5.01 / 9.28 s |
| Grant one 40-digit ID | 3.51 / 6.04 s | 6.54 / 12.82 s |
| Revoke the ID | 4.72 / 8.41 s | 7.96 / 15.94 s |

Both runs used x86-64 Linux, local NVMe storage with ext4, an Intel Core i5-12450H, GMP 6.3.0 and FLINT 3.0.1. Timings run from subprocess launch through completion and include the normal file and directory syncs on changed writes. With warm input, kernel-accounted storage reads were at most 20 KiB per operation. After `POSIX_FADV_DONTNEED`, each call incurred about 2 GiB of storage reads. Hardware caches remained under device control. The largest sampled post-exec peak resident memory across these twelve operations was **6,644 KiB**. The memory cap was enforced independently of sampling. Kernel page cache and disk space for the replacement file are additional resources. [Raw results from both runs](validation.json)

## Reproduce

Run these commands from the repository root on Linux with Python 3.11 or later.

```sh
cargo build --release --locked
python3 tools/verify_scale.py --checkpoint --size-mib 2048 --output /tmp/prime-label-results.json
```

Allow slightly more than 4 GiB of scratch space on the filesystem to measure. Select it with `--scratch DIRECTORY`. The harness regenerates the historical member stream from seed 20261128, checks its known digest, then passes the IDs through the public CLI's primality proofs. For a quick smoke test omit `--checkpoint` and use `--size-mib 1`. [Verification harness](../tools/verify_scale.py)
