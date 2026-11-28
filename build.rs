fn main() {
    println!("cargo:rerun-if-changed=src/prime_check.c");
    let flint = pkg_config::Config::new()
        .cargo_metadata(false)
        .probe("flint")
        .expect("install the FLINT development package and pkg-config");
    cc::Build::new()
        .file("src/prime_check.c")
        .includes(flint.include_paths)
        .compile("prime_check");
    for path in flint.link_paths {
        println!("cargo:rustc-link-search=native={}", path.display());
    }
    // The Ubuntu FLINT 3.0.1 .pc file lists its dependencies but omits FLINT.
    println!("cargo:rustc-link-lib=flint");
}
