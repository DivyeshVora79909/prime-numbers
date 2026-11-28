# Prime labels

## Decision record

An entity has an externally assigned prime ID and a label representing the set of IDs it admits. The application owns entity names, identity allocation and the choice of which resource to query. This library owns one exact integer per label. An empty set has product 1. Adding an ID introduces one prime factor, and removing it eliminates that factor. Unique factorization makes membership equivalent to a zero remainder.

The previous project optimized signed recursive powers and demonstrated exact arithmetic. At 30,000 random 40-digit primes its H2N1 label needed 15,686,193 bytes, already within 91 bytes of that format's optimum. Ordinary binary needed 494,227 bytes. Packing nonadjacent signed digits also needed 494,227 payload bytes. In the controlled 64 KiB block experiment, binary had simpler arithmetic and file updates were faster on that checkpoint. Sparse towers were a genuinely different workload. The new product is specifically for dense prime-set labels. [Historical evidence and source](archive/README.md)

For a product `N`, the payload uses `ceil((floor(log2(N)) + 1) / 8)` bytes. The logarithm of a product is the sum of the logarithms of its factors, so each new 40-digit prime adds about 130–133 bits to that sum. Each membership check reads the resulting file. This makes both storage growth and query cost predictable as a set grows.

## Smallest useful public contract

* A label is a nonempty file containing a positive, shortest unsigned little-endian magnitude. The most significant byte is nonzero. The empty set is the single byte `01`. Equal sets produce byte-identical files regardless of update history. A file produced outside the API must also satisfy the squarefree product invariant.
* The Rust library and a thin `prime-label` CLI expose initialization, membership, remainder, and one declarative batch update. Callers describe grants and revocations. A repeated grant or absent revocation is a no-op, identical requests deduplicate, and conflicting directions in a batch are rejected before publication.
* Individual entity IDs are decimal primes, including 2. Use FLINT's primality proof at the ID boundary. Generate a new 40-digit ID from OS randomness or import an existing prime. Uniqueness of entity assignments belongs to the caller's registry. The general remainder operation accepts any positive modulus.
* Store only the product. Readers open one immutable snapshot. Writers serialize on a stable per-label lock, use a temporary file in the same directory, flush its contents, rename it atomically, then sync the directory. Retryable set semantics cover uncertainty after a successful rename. Initialization must never overwrite a label.

## Implementation plan

1. Freeze the exact prior tree, the uncommitted benchmark work, and the supplied research transcript in a compressed, verified archive. Keep a short index and restoration command.
2. Use one Rust crate. A file arithmetic module owns block boundaries and calls Rug/GMP, a prime ID type owns validation, and the library owns batch and publication semantics. A small CLI delegates to these functions. Use maintained libraries for integer arithmetic, primality and temporary files. No codec registry, symbolic grammar or optimizer belongs in the active path.
3. Given the distinct requested grant product A and removal product D, revoke gcd(N,D) and grant A/gcd(N,A). Compute these intersections from streamed remainders. Build large batch products with balanced multiplication. Size working memory by the block and batch, independent of the current label's length for a fixed request.
4. Scan from the high end for remainder and exact division. Scan from the low end for multiplication. Mixed updates may use the unpublished quotient as the in-place multiplication source. An error before publication leaves the prior generation available. A single fixed-width prime costs linear work in the number of stored bytes. Bulk arithmetic also depends on total operand width.
5. Expose build, examples and verification from the repository root. Retire old source trees, package entry points, generated artifacts and duplicate documentation after checking the archive.

## Acceptance criteria

* Test mathematical results against an independent implementation, including 2, random 40-digit IDs, absent/present duplicates, conflicting edits, malformed labels, chunk boundaries and mixed batches. A failed update preserves the previous bytes. Concurrent writers cannot lose accepted grants.
* Rebuild the known 30,000-prime label and match its packed SHA-256. Verify the independent outsider remainder, remove and restore the member, and compare the full restored file.
* Exercise public commands against a real, densely written input of at least 1 GiB under a smaller process memory limit. Verify publication and an exact round trip by streaming hashes. Distinguish synthetic scale from prime provenance, buffered cache effects from storage durability, and measurement from any claim of global optimality.
* Verify a clean build, tests, formatting, lints, runnable CLI help, and an independently restorable archive.
