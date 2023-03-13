use std::{fs::canonicalize, path::Path};

use assert_cmd::Command;
use assert_fs::fixture::NamedTempFile;
use predicates::prelude::*;

use rust_htslib::bam::{self, Read, Records};

macro_rules! canon_arg {
    (flag $flag:expr)       => {[$flag]};
    ($flag:expr, $val:expr) => {[$flag, canon_arg!(val $val)] };
    (val $val:expr)         => {
        canonicalize($val)
            .expect(&format!("Failed to canonicalize {}", $val))
            .to_str()
            .expect("Invalid path")
    };
}

fn rust_htslib_read_back(output: &Path) -> bam::Reader {
    bam::Reader::from_path(&output).expect("Failed to open output back.")
}

fn output_is_masked(output: &Path) ->  bool {
    rust_htslib_read_back(output).records().any(|rec| rec.expect("Invalid Record").seq().as_bytes().contains(&b'N'))
}

#[test]
fn basic_command_stdout() {
    let cmd = Command::cargo_bin("pmd-mask").expect("Invalid")
    .args(canon_arg!("--reference", "tests/test-data/reference/hs37d5-MTonly/hs37d5-MTonly.fa.gz"))
    .args(canon_arg!("--bam", "tests/test-data/bam/dummy-MTonly/dummy-MTonly-1000.bam"))
    .args(canon_arg!("--misincorporation", "tests/test-data/bam/dummy-MTonly/misincorporation.txt"))
    .assert();

    println!("{cmd}");

    cmd.success()
    .code(0)
    .stderr(predicate::str::is_empty()) // No Logging warnings, errors, etc.
    .stdout(predicate::str::starts_with("@HD"));


}

#[test]
fn basic_command_samfile_output() {

    let fixture_bam = NamedTempFile::new("output.bam").expect("Failed to create fixture for output bam");
    println!("Output bam fixture: {}", fixture_bam.display());

    let cmd = Command::cargo_bin("pmd-mask").expect("Invalid")
    .args(canon_arg!("--reference", "tests/test-data/reference/hs37d5-MTonly/hs37d5-MTonly.fa.gz"))
    .args(canon_arg!("--bam", "tests/test-data/bam/dummy-MTonly/dummy-MTonly-1000.bam"))
    .args(canon_arg!("--misincorporation", "tests/test-data/bam/dummy-MTonly/misincorporation.txt"))
    .args(["--output", fixture_bam.to_str().expect("Non UTF8 character in fixture")])
    .assert();

    println!("{cmd}");

    // ---- Check if the run was successful
    cmd.success()
        .code(0)
        .stderr(predicate::str::is_empty())  // No warnings, no errors, no info
        .stdout(predicate::str::is_empty()); // Since we're writing to file, stdout should be empty.

    // ---- rust_htslib should be able to read the entire contents. We also assume some masking applied (since the default is 0.01)
    assert!(output_is_masked(&fixture_bam));
    

    fixture_bam.close().expect("Failed to delete fixture");
}

#[test]
fn complex_command() {
    let fixture_bam = NamedTempFile::new("output.bam").expect("Failed to create fixture for output bam");
    println!("Output bam fixture: {}", fixture_bam.display());

    let cmd = Command::cargo_bin("pmd-mask").expect("Invalid")
    .args(canon_arg!("--reference", "tests/test-data/reference/hs37d5-MTonly/hs37d5-MTonly.fa.gz"))
    .args(canon_arg!("--bam", "tests/test-data/bam/dummy-MTonly/dummy-MTonly-1000.bam"))
    .args(canon_arg!("--misincorporation", "tests/test-data/bam/dummy-MTonly/misincorporation.txt"))
    .args(["--output", fixture_bam.to_str().expect("Non UTF8 character in fixture")])
    .args(["--threshold", "1.0"])
    .args(["--output-fmt", "BAM"])
    .args(["-vvv"])
    .assert();

    println!("{cmd}");

    cmd.success()
        .code(0)
        .stderr(predicate::str::contains("INFO"))  // Contains INFO level logging information
        .stderr(predicate::str::contains("DEBUG")) // Contains DEBUG level
        .stderr(predicate::str::contains("TRACE")) // Contains TRACE level 
        .stdout(predicate::str::is_empty()); // Since we're writing to file, stdout should be empty.

    assert!(!output_is_masked(&fixture_bam));

    fixture_bam.close().expect("Failed to delete fixture");
}
