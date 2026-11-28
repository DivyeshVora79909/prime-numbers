use prime_label::{Changes, Integer, PrimeId, apply, contains, init, remainder};
use std::fs::{self, File};
use std::io::{Read, Write};
use std::os::unix::fs::{MetadataExt, PermissionsExt, symlink};
use std::process::{Command, Stdio};
use tempfile::tempdir;

fn id(n: u32) -> PrimeId {
    n.to_string().parse().unwrap()
}

#[test]
fn batches_are_exact_sets_and_leave_readers_on_their_snapshot() {
    let work = tempdir().unwrap();
    let path = work.path().join("access");
    init(&path).unwrap();
    let mut old = File::open(&path).unwrap();
    fs::set_permissions(&path, fs::Permissions::from_mode(0o640)).unwrap();

    let mut edits = Changes::default();
    for p in [2, 3, 5, 7, 7] {
        edits.grant(id(p)).unwrap();
    }
    assert!(apply(&path, edits).unwrap());
    assert_eq!(fs::read(&path).unwrap(), [210]);
    let mut snapshot = Vec::new();
    old.read_to_end(&mut snapshot).unwrap();
    assert_eq!(snapshot, [1]);
    assert_eq!(path.metadata().unwrap().permissions().mode() & 0o777, 0o640);

    let changed = || {
        let mut batch = Changes::default();
        for p in [2, 5, 19] {
            batch.revoke(id(p)).unwrap();
        }
        for p in [3, 11, 13] {
            batch.grant(id(p)).unwrap();
        }
        batch
    };
    assert!(apply(&path, changed()).unwrap());
    let inode = path.metadata().unwrap().ino();
    assert!(!apply(&path, changed()).unwrap());
    assert_eq!(path.metadata().unwrap().ino(), inode);
    assert_eq!(fs::read(&path).unwrap(), (3u16 * 7 * 11 * 13).to_le_bytes());
    assert_eq!(remainder(&path, &Integer::from(28)).unwrap(), 3003 % 28);
    assert_eq!(remainder(&path, &Integer::from(1)).unwrap(), 0);
    assert!(contains(&path, &id(7)).unwrap());
    assert!(!contains(&path, &id(2)).unwrap());
    let mut all = Changes::default();
    for p in [3, 7, 11, 13] {
        all.revoke(id(p)).unwrap();
    }
    assert!(apply(&path, all).unwrap());
    assert_eq!(fs::read(&path).unwrap(), [1]);
}

#[test]
fn invalid_inputs_do_not_replace_existing_labels() {
    let work = tempdir().unwrap();
    let path = work.path().join("label");
    init(&path).unwrap();
    assert!(init(&path).is_err());
    assert!(init(work.path().join("label.lock")).is_err());
    let alias = work.path().join("link");
    symlink(&path, &alias).unwrap();
    let mut add = Changes::default();
    add.grant(id(2)).unwrap();
    assert!(apply(&alias, add).is_err());
    assert!(alias.is_symlink());
    assert!(remainder(&path, &Integer::from(0)).is_err());
    assert_eq!(fs::read(&path).unwrap(), [1]);
    for bad in [&[][..], &[0], &[1, 0]] {
        fs::write(&path, bad).unwrap();
        assert!(apply(&path, Changes::default()).is_err());
        assert!(contains(&path, &id(7)).is_err());
        assert_eq!(fs::read(&path).unwrap(), bad);
    }
}

#[test]
fn competing_cli_processes_preserve_every_grant() {
    let work = tempdir().unwrap();
    let path = work.path().join("label");
    init(&path).unwrap();
    let guard = File::open(work.path().join("label.lock")).unwrap();
    guard.lock().unwrap();
    let primes = [2u32, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    let writers: Vec<_> = primes
        .iter()
        .map(|p| {
            Command::new(env!("CARGO_BIN_EXE_prime-label"))
                .arg("apply")
                .arg(&path)
                .arg(format!("+{p}"))
                .stdout(Stdio::null())
                .spawn()
                .unwrap()
        })
        .collect();
    guard.unlock().unwrap();
    for mut child in writers {
        assert!(child.wait().unwrap().success());
    }
    let n: u64 = primes.iter().map(|&p| u64::from(p)).product();
    let bytes = fs::read(&path).unwrap();
    assert_eq!(Integer::from_digits(&bytes, rug::integer::Order::Lsf), n);
}

#[test]
fn textual_batches_parse_entirely_before_mutation() {
    let work = tempdir().unwrap();
    let path = work.path().join("label");
    init(&path).unwrap();
    for (text, success, expected) in [
        ("+11 +13\n+11\n", true, 143),
        ("-11 +17\n+4\n", false, 143),
        ("+17 -13 -17", false, 143),
        ("-13 +2", true, 22),
    ] {
        let mut cli = Command::new(env!("CARGO_BIN_EXE_prime-label"));
        let mut child = cli
            .arg("apply")
            .arg(&path)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        child
            .stdin
            .take()
            .unwrap()
            .write_all(text.as_bytes())
            .unwrap();
        assert_eq!(child.wait_with_output().unwrap().status.success(), success);
        assert_eq!(fs::read(&path).unwrap(), [expected]);
    }
}
