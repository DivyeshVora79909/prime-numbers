"""Exercise the installed CLI on Linux with a 96 MiB address-space limit."""

import argparse
import hashlib
import json
import math
import os
import platform
import resource
import select
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from random import Random

MIB = 1 << 20
BLOCK = 64 << 10
MEMBER = 3589358312608496142804550037657486194871
OUTSIDER = 1540147111665682694523005535063805113443
IDS_HASH = "39c157d6c06b8794421b7ca1bf884d2faebc909de1d23637f968918aba82e9f5"
LABEL_HASH = "b4a6c3bb41f08ef79c63ab6cbc0dcaa65db2405c6ace6819859366ec80021c5b"
QUOTIENT_HASH = "e9802545f34d426a11578449eb268a9e8031200d70cff958cf91cbc9d4251991"
OUTSIDER_REMAINDER = 833120395990529729763735569905656999601


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def limit_worker():
    resource.setrlimit(resource.RLIMIT_AS, (96 * MIB, 96 * MIB))


def run(binary, label, operation, *args, stdin=b""):
    """Time through exit. Linux wait4 reports I/O in units of 512 bytes."""
    with tempfile.TemporaryFile() as src, tempfile.TemporaryFile() as out:
        src.write(stdin)
        src.seek(0)
        start = time.perf_counter()
        with subprocess.Popen([binary, operation, label, *map(str, args)], stdin=src,
                              stdout=out, stderr=out, preexec_fn=limit_worker) as child:
            proc = Path(f"/proc/{child.pid}")
            hwm = 0
            pidfd = os.pidfd_open(child.pid)
            try:
                while not select.select([pidfd], [], [], 0.02)[0]:
                    try:
                        if (proc / "exe").readlink() == binary:
                            for line in (proc / "status").read_text().splitlines():
                                if line.startswith("VmHWM:"):
                                    hwm = max(hwm, int(line.split()[1]))
                    except (FileNotFoundError, ProcessLookupError):
                        pass  # Exited between pidfd polling and /proc access.
                _, status, usage = os.wait4(child.pid, 0)
                child.returncode = os.waitstatus_to_exitcode(status)
            finally:
                os.close(pidfd)
        seconds = time.perf_counter() - start
        out.seek(0)
        output = out.read().decode().strip()
        if child.returncode:
            raise RuntimeError(f"{operation} exited {child.returncode}: {output}")
    return output, {"seconds": round(seconds, 6), "peak_rss_kib_including_launcher": usage.ru_maxrss,
                    "observed_post_exec_hwm_kib": hwm or None,
                    "storage_read_bytes": usage.ru_inblock * 512, "storage_write_bytes": usage.ru_oublock * 512}


def checkpoint_bytes():
    """Replay the archived draw order; the pinned hash identifies the entire stream."""
    rng = Random(20261128)
    small_primes = math.prod((3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47))
    seen, records = set(), bytearray()
    while len(records) < 30_000 * 17:
        p = rng.randrange(10**39 + 1, 10**40, 2)
        if p in seen or math.gcd(p, small_primes) != 1:
            continue
        # Cheap candidate filter. The release CLI proves every ID with FLINT.
        if pow(2, p - 1, p) != 1 or pow(3, p - 1, p) != 1:
            continue
        index = len(seen)
        seen.add(p)
        if index not in (0, 3) and not 12292 <= index < 12324:
            records.extend(p.to_bytes(17, "big"))
    return records


def actual_primes(binary, label, raw):
    if len(raw) != 30_000 * 17 or hashlib.sha256(raw).hexdigest() != IDS_HASH:
        raise ValueError("30,000-prime candidate stream differs from the archived checkpoint")
    ids = [int.from_bytes(raw[i:i + 17], "big") for i in range(0, len(raw), 17)]
    assert len(set(ids)) == len(ids) and ids[-1] == MEMBER
    run(binary, label, "init")
    requests = "".join(f"+{p}\n" for p in ids).encode()
    result, build = run(binary, label, "apply", stdin=requests)
    assert result == "changed"
    assert label.stat().st_size == 494_227 and digest(label) == LABEL_HASH
    assert run(binary, label, "contains", MEMBER)[0] == "true"
    assert run(binary, label, "contains", OUTSIDER)[0] == "false"
    assert int(run(binary, label, "mod", OUTSIDER)[0]) == OUTSIDER_REMAINDER
    run(binary, label, "apply", f"-{MEMBER}", f"-{MEMBER}", f"-{OUTSIDER}")
    assert digest(label) == QUOTIENT_HASH and run(binary, label, "contains", MEMBER)[0] == "false"
    run(binary, label, "apply", f"+{MEMBER}", f"+{MEMBER}")
    assert run(binary, label, "apply", f"+{MEMBER}", f"-{OUTSIDER}")[0] == "unchanged"
    assert digest(label) == LABEL_HASH
    return {"ids": len(ids), "bytes": label.stat().st_size, "sha256": LABEL_HASH, "bulk_build": build}


def synthetic_file(label, size):
    """Generate dense arithmetic input; primality of its other factors is unspecified."""
    granted_hash = hashlib.sha256()
    carry = 0
    radix_mask = (1 << (8 * BLOCK)) - 1
    with label.open("wb") as output:
        for offset in range(0, size, BLOCK):
            data = bytearray(hashlib.shake_256(f"prime-label/scale/v1/{size}/{offset}".encode()).digest(BLOCK))
            if offset == 0:
                data[0] |= 1
            if offset + BLOCK == size:
                data[-1] |= 128
            output.write(data)
            value = int.from_bytes(data, "little") * MEMBER + carry
            granted_hash.update((value & radix_mask).to_bytes(BLOCK, "little"))
            carry = value >> (8 * BLOCK)
        output.flush()
        os.fsync(output.fileno())
    granted_hash.update(carry.to_bytes((carry.bit_length() + 7) // 8, "little"))
    return digest(label), granted_hash.hexdigest()


def remainder(path, modulus):
    with path.open("rb") as stream:
        end = stream.seek(0, os.SEEK_END)
        r = 0
        while end:
            n = min(BLOCK, end)
            end -= n
            stream.seek(end)
            r = ((r << (8 * n)) + int.from_bytes(stream.read(n), "little")) % modulus
        return r


def prepare_cache(label, mode):
    with label.open("rb", buffering=0) as stream:
        if mode == "warm":
            buf = bytearray(MIB)
            while stream.readinto(buf):
                pass
        else:
            os.posix_fadvise(stream.fileno(), 0, 0, os.POSIX_FADV_DONTNEED)


def dense_scale(binary, label, size):
    original, granted = synthetic_file(label, size)
    assert remainder(label, MEMBER) != 0
    expected = remainder(label, OUTSIDER)
    assert int(run(binary, label, "mod", OUTSIDER)[0]) == expected
    cases = {}
    for mode in ("warm", "eviction_requested"):
        operations = {}
        for name, operation, arg in [("grant", "apply", f"+{MEMBER}"), ("member", "contains", MEMBER),
                                     ("revoke", "apply", f"-{MEMBER}")]:
            prepare_cache(label, mode)
            answer, operations[name] = run(binary, label, operation, arg)
            if name == "member":
                assert answer == "true"
            elif name == "grant":
                assert digest(label) == granted
                assert int(run(binary, label, "mod", OUTSIDER)[0]) == expected * MEMBER % OUTSIDER
                if mode == "warm":
                    assert run(binary, label, "apply", f"+{MEMBER}")[0] == "unchanged"
                    assert digest(label) == granted
            else:
                assert digest(label) == original and run(binary, label, "contains", MEMBER)[0] == "false"
            print(f"{mode} {name}: {operations[name]['seconds']:.3f}s", flush=True)
        cases[mode] = operations
    return {"fixture": "SHAKE-256 dense arithmetic cofactor; other factors uncertified",
            "cofactor_bytes": size, "cofactor_sha256": original, "granted_sha256": granted, "cases": cases}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--binary", type=Path, default=Path("target/release/prime-label"))
    checkpoint = ap.add_mutually_exclusive_group()
    checkpoint.add_argument("--checkpoint", action="store_true", help="regenerate the seeded 30,000-prime checkpoint")
    checkpoint.add_argument("--members", type=Path, help="use a saved checkpoint of 17-byte big-endian IDs")
    ap.add_argument("--size-mib", type=int, default=1, help="actual synthetic file size (default: smoke test)")
    ap.add_argument("--scratch", type=Path, help="filesystem used for temporary label and replacement")
    ap.add_argument("--output", type=Path, help="also save the JSON report to a file")
    args = ap.parse_args()
    if not __debug__ or args.size_mib < 1:
        ap.error("enable Python assertions and select a positive MiB count")
    binary = args.binary.resolve(strict=True)
    size = args.size_mib * MIB
    with tempfile.TemporaryDirectory(prefix="prime-label-check-", dir=args.scratch) as temp:
        directory = Path(temp)
        if shutil.disk_usage(directory).free < 2 * size + 64 * MIB:
            ap.error("insufficient scratch space for the label and its atomic replacement")
        report = {"platform": platform.platform(), "worker_address_space_limit_mib": 96,
                  "executable_sha256": digest(binary)}
        if args.members or args.checkpoint:
            raw = args.members.read_bytes() if args.members else checkpoint_bytes()
            report["actual_primes"] = actual_primes(binary, directory / "actual", raw)
        report["synthetic"] = dense_scale(binary, directory / "synthetic", size)
        if digest(binary) != report["executable_sha256"]:
            raise RuntimeError("release executable changed during validation")
    text = json.dumps(report, indent=2)
    print(text)
    if args.output:
        args.output.write_text(text + "\n")


if __name__ == "__main__":
    main()
