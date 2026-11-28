#![doc = include_str!("../README.md")]

mod arithmetic;
mod prime;

pub use prime::PrimeId;
pub use rug::Integer;

use std::collections::{BTreeSet, VecDeque};
use std::fs::{File, OpenOptions};
use std::io::{self, Write};
use std::os::unix::fs::OpenOptionsExt;
use std::path::Path;
use tempfile::NamedTempFile;

/// A declarative, retryable set update. Repetitions in one direction deduplicate.
#[derive(Default)]
pub struct Changes {
    grants: BTreeSet<PrimeId>,
    revocations: BTreeSet<PrimeId>,
}

impl Changes {
    pub fn grant(&mut self, id: PrimeId) -> io::Result<()> {
        insert(&mut self.grants, &self.revocations, id)
    }

    pub fn revoke(&mut self, id: PrimeId) -> io::Result<()> {
        insert(&mut self.revocations, &self.grants, id)
    }
}

fn insert(
    into: &mut BTreeSet<PrimeId>,
    opposite: &BTreeSet<PrimeId>,
    id: PrimeId,
) -> io::Result<()> {
    if opposite.contains(&id) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("conflicting edits for {id}"),
        ));
    }
    into.insert(id);
    Ok(())
}

/// Combine similarly sized inputs first. A batch of k IDs uses k−1 products.
fn product(ids: BTreeSet<PrimeId>) -> Integer {
    let mut pending: VecDeque<_> = ids.into_iter().map(PrimeId::into_integer).collect();
    while pending.len() > 1 {
        let mut a = pending.pop_front().unwrap();
        a *= pending.pop_front().unwrap();
        pending.push_back(a);
    }
    pending.pop_front().unwrap_or_else(|| Integer::from(1))
}

fn parent(path: &Path) -> &Path {
    path.parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."))
}

fn lock(path: &Path) -> io::Result<File> {
    let name = path
        .file_name()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "expected a label filename"))?;
    if name.as_encoded_bytes().ends_with(b".lock") {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "the .lock suffix is reserved",
        ));
    }
    let mut sidecar = name.to_os_string();
    sidecar.push(".lock");
    let path = parent(path).join(sidecar);
    let file = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(false)
        .mode(0o600)
        .custom_flags(libc::O_NOFOLLOW | libc::O_NONBLOCK)
        .open(path)?;
    if !file.metadata()?.is_file() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "expected a regular lock file",
        ));
    }
    file.lock()?;
    Ok(file)
}

fn snapshot(path: &Path) -> io::Result<File> {
    OpenOptions::new()
        .read(true)
        .custom_flags(libc::O_NOFOLLOW | libc::O_NONBLOCK)
        .open(path)
}

fn publish(tmp: NamedTempFile, path: &Path, create: bool) -> io::Result<()> {
    tmp.as_file().sync_all()?;
    if create {
        tmp.persist_noclobber(path)?;
    } else {
        tmp.persist(path)?;
    }
    File::open(parent(path))?.sync_all()
}

/// Create an empty set. An existing file is preserved.
pub fn init(path: impl AsRef<Path>) -> io::Result<()> {
    let path = path.as_ref();
    let _guard = lock(path)?;
    let mut tmp = NamedTempFile::new_in(parent(path))?;
    tmp.write_all(&[1])?;
    publish(tmp, path, true)
}

/// Read one snapshot and return its remainder for any positive modulus.
pub fn remainder(path: impl AsRef<Path>, modulus: &Integer) -> io::Result<Integer> {
    arithmetic::remainder(&snapshot(path.as_ref())?, modulus)
}

pub fn contains(path: impl AsRef<Path>, id: &PrimeId) -> io::Result<bool> {
    Ok(remainder(path, id.value())? == 0)
}

/// Atomically apply all grants and revocations. True means the set changed.
pub fn apply(path: impl AsRef<Path>, changes: Changes) -> io::Result<bool> {
    let path = path.as_ref();
    let mut add = product(changes.grants);
    let mut remove = product(changes.revocations);
    let _guard = lock(path)?;
    let input = snapshot(path)?;
    // The requested products are coprime. One remainder gives both intersections.
    let combined = Integer::from(&add * &remove);
    let residue = arithmetic::remainder(&input, &combined)?;
    let already_present = residue.clone().gcd(&add);
    add.div_exact_mut(&already_present);
    remove.gcd_mut(&residue);
    if add == 1 && remove == 1 {
        // A retry can observe a published generation whose directory sync failed.
        File::open(parent(path))?.sync_all()?;
        return Ok(false);
    }

    let tmp = NamedTempFile::new_in(parent(path))?;
    let out = tmp.as_file();
    if remove != 1 {
        arithmetic::divide(&input, &remove, out)?;
        if add != 1 {
            arithmetic::multiply(out, &add, out)?;
        }
    } else {
        arithmetic::multiply(&input, &add, out)?;
    }
    out.set_permissions(input.metadata()?.permissions())?;
    publish(tmp, path, false)?;
    Ok(true)
}
