# Contributing

## Running benchmarks and regression testing
Listing benchmarrks: 
```Bash
cargo bench -- --list 2>&1 | grep "bench$"
```

Running a quick bench using all benchmark files (3s warmup time, followed by 5s measurement time)
```Bash
cargo bench
```

An HTML report should be available in the `target/criterion` subdirectory
```Bash
firefox ./target/criterion/report/index.html
```

Running a specific bench: 
```Bash
cargo bench --bench apply_pmd_mask -- --warm-up-time 10 --measurement-time 60
```

## Getting code coverage metrics for the pmd-mask codebase
Install [llvm-cov](https://github.com/taiki-e/cargo-llvm-cov) 
```Bash
cargo install llvm-cov
```

Run `llvm-cov`
```Bash
# Using shell output. either stdout or pager
cargo llvm-cov --workspace --all 
# or 
cargo llvm-cov --worskpace --all --text | less -R

# Html report
cargo llvm-cov --workspace --all --open

# lcov format
cargo llvm-cov --workspace --all --lcov > lcov.info
```

## Building releases for cross-compilation
the use of [`cross-rs`](https://github.com/cross-rs/cross) is here recommended for seamless cross-compilation of binaries. Note that `cross-rs` will generally require administrative privileges, as it makes use of `docker` for containerization. `podman` can be used as an alternative to docker in cases where this is not possible. See the detailled `cross-rs` [dependencies](https://github.com/cross-rs/cross?tab=readme-ov-file#dependencies), [installation instructions](https://github.com/cross-rs/cross?tab=readme-ov-file#installation), and [wiki](https://github.com/cross-rs/cross/wiki) for additional information. In essence, and assuming default configuration, `cross-rs` can be installed using the following cargo command:
```bash
cargo install cross --git https://github.com/cross-rs/cross
```

Next, you'll need to ensure the `dockerd` and `containerd` are properly running on your system before continuing. (Note this only needs to be checked once per boot):
```bash
sudo systemctl start docker
```

Cross can then be executed to build, and test cross-compiled binaries:
```bash
TARGET_ARCH="x86_64-unknown-linux-gnu"
cross test --target ${TARGET_ARCH} && cross build --release --target ${TARGET_ARCH}
```

## Checking the minimum supported rust version (MSRV) of the software:

[`cargo-msrv`](https://github.com/foresterre/cargo-msrv?tab=readme-ov-file) is here recommended to quickly check the MSRV of `pmd-mask` for source compilation. If you wish to submit a pull-request to `pmd-mask`, please check and report any changes in MSRV, resulting from your changes.

- **<ins>Installation</ins>**: `cargo install cargo-msrv --locked`
- **<ins>Usage</ins>**: `cargo msrv find`

