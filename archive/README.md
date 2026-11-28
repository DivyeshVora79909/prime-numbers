# Before packed-binary labels

The archive captures the working project immediately before the redesign. Its inventory consists of the tracked files at `2b2b57b`, the six additional packed-scale benchmark files from that checkout, and the original user-supplied research conversation. The snapshot preserves source paths, tests, benchmark data and the storage-growth figure. The bundle uses sorted relative paths, fixed ownership and timestamps, and XZ compression. Verification compares every extracted file with its source.

## Research in three steps

1. The first demonstration answered divisibility queries over retained prime factors. Recursive signed base-2 arithmetic then merged those factors into one integer and supported exact division and modular queries on its saved H2N1 circuit.
2. At 30,000 random 40-digit primes, the H2N1 label occupied 15,686,193 bytes. Its format-specific lower bound was 15,686,102 bytes. The same positive integer needed 494,227 bytes in ordinary packed binary. The earlier Rust kernels accelerated arithmetic and serialization while preserving the full H2N1 output.
3. A packed signed-digit stream achieved the same 494,227-byte payload. With identical GMP arithmetic, the file-backed comparison measured binary / packed NAF medians of 1.88 / 2.21 ms for membership, 4.57 / 5.04 ms for removal and 5.20 / 6.15 ms for insertion on this checkpoint. Both restored the original bytes. Inputs were warmed, and writes were buffered. The recorded synthetic run reached 128 MiB under a 96 MiB address-space cap. Sparse towers and geometric expressions remained strong cases for symbolic storage. Dense access labels favored ordinary binary.

## Reading the snapshot

| Path after extraction | Material |
| --- | --- |
| `README.md`, `docs/` | Historical API, algorithm and proof details, 40-digit prime sampler, checkpoint reproduction, graph and measurements |
| `src/symbolic_primes/`, `native/`, `tests/` | Python engine, Rust extension and signed-arithmetic regression tests |
| `benchmarks/results/` | Raw benchmark outputs, 17,710-prefix storage sweep and packed-scale results |
| `benchmarks/packed_scale/` | The uncommitted matched-kernel binary/NAF comparison and its verification harness |
| `context/Primality Testing After RH.txt` | Verbatim chat supplied by the user. It records research questions and proposals. The benchmark files document the results verified in the project. |

## Restore and verify

The bundle is 797,092 bytes and contains 45 tracked files, six additional benchmark files, one transcript and a SHA-256 manifest. The 52 original files total 2,603,999 bytes. The transcript's SHA-256 is `f4e3fe49bc5a0e59bd359745b5e8837c406218a2d20ecc383c39cad05fea9cd3`.

From the project root, restore to a fresh directory:

```sh
printf '%s  %s\n' fb70372dfcaaf5ac0c161716e7aa13e3ffc0f036f2e4d5a6f680dd83ee1d4c56 archive/before-packed-binary.tar.xz | sha256sum --check
snapshot=$(mktemp -d)
tar -xJf archive/before-packed-binary.tar.xz -C "$snapshot"
(cd "$snapshot" && sha256sum --check SNAPSHOT.sha256)
```

The original tracked tree is also available at Git commit `2b2b57b1fcbd67202471db8315ef7de4da157c9f`. The snapshot includes the six later benchmark files at their original paths. The archived reports include recipes for regenerating the large checkpoint files.

To reproduce the archive hash from the extracted files with GNU tar and XZ:

```sh
(cd "$snapshot" && find . -type f -printf '%P\0' | LC_ALL=C sort -z | tar --null -T - --format=ustar --owner=0 --group=0 --numeric-owner --mtime=@0 --mode=0644 -cf -) | xz -9e --threads=1 --check=crc64 -c | sha256sum
```
