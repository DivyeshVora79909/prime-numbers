# Prime labels

An entity has a prime ID. A resource's access label is the product of the prime IDs it admits. To check membership, divide the label by the presented ID and inspect the remainder. The empty set is 1.

`prime-label` stores the product as a single packed binary file. Rust streams it in 64 KiB blocks, GMP computes the arithmetic, and FLINT proves primality when an ID is imported or generated. A batch of grants and revocations publishes one new generation of the file.

## Build and try it

On Ubuntu 24.04 with rustup installed, build the executable and put it on your path.

```sh
sudo apt-get install build-essential pkg-config libgmp-dev libflint-dev
cargo build --release --locked
export PATH="$PWD/target/release:$PATH"
```

Run this example in a new directory.

```sh
alice=$(prime-label new-id)
prime-label init resource.label
prime-label apply resource.label "+$alice" +101
prime-label contains resource.label "$alice"    # true
prime-label apply resource.label "-$alice" +103
prime-label contains resource.label "$alice"    # false
prime-label mod resource.label 103             # 0
```

`new-id` returns a random proven 40-digit prime. Keep the ID with the entity record and enforce uniqueness in your registry. The membership call prints `true` or `false`. A command failure exits with a nonzero status.

For bulk edits, pipe whitespace-separated signed IDs to `prime-label apply resource.label`. Additions of present IDs and removals of absent IDs leave the set as requested. Repetitions deduplicate. Asking for both directions on one ID rejects the entire batch. Successful edits print `changed` or `unchanged`.

## Rust API

```rust
use prime_label::{Changes, PrimeId, apply, contains, init};

fn main() -> std::io::Result<()> {
    let dir = tempfile::tempdir()?;
    let label = dir.path().join("resource.label");
    let alice: PrimeId = "3589358312608496142804550037657486194871".parse()?;
    init(&label)?;
    let mut changes = Changes::default();
    changes.grant(alice.clone())?;
    assert!(apply(&label, changes)?);
    assert!(contains(&label, &alice)?);
    Ok(())
}
```

## Storage and operating cost

A label contains the shortest positive unsigned little-endian magnitude. For example, granting 3 and 5 produces one byte, `0f`. A zero most significant byte or an empty file is invalid. Imported files must encode squarefree products for set semantics. Labels created and edited by this library preserve that invariant.

Use an application-owned directory on a local Linux filesystem. One stable sibling `LABEL.lock` serializes writers. The suffix is reserved. Readers open one snapshot. Changed updates require room for the replacement file, preserve permission bits, and sync both the file and its directory. Keep the sidecar in place across updates. A directory-sync error after rename can leave the new generation visible. Repeating the same batch is safe.

For a single fixed-width ID, membership and updates take linear time in the file size and bounded working memory. A changed update first checks membership and then streams the arithmetic. A mixed grant/revoke can reuse the temporary quotient in place. Bulk memory also depends on the IDs supplied and the size of their combined products.

## Verification

```sh
cargo test --locked
cargo clippy --locked --all-targets -- -D warnings
python3 tools/verify_scale.py --size-mib 1
```

The 30,000-prime checkpoint and the 2 GiB disk test have a self-contained reproduction command. [Measured results](docs/validation.md) · [Design and decisions](DESIGN.md) · [Historical archive and restore instructions](archive/README.md)
