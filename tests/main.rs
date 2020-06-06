use assert_cmd::prelude::*;
use std::process::Command;
use std::process::ExitStatus;

fn sc_jaccard() -> std::process::Command { 
    Command::cargo_bin("scJaccard").unwrap()
}

#[test]
fn help() {
    sc_jaccard().arg("-h").assert().success();
}

#[test]
fn test_success_run() { 
    let out = sc_jaccard()
        .args(&["-i", "example_data/test.mtx", "--atac", "example_data/target.bed.gz", "--barcodes", "example_data/test.tsv", "--nchr", "1", "--loglevel", "warn", "--bed", "example_data/test.bed"])
        .assert();

    let output = out
        .get_output();
    let known_out = "barcode2 0.14286\nbarcode1 0.33333\n";

    println!("{:?}", output);
    println!("{:?}", out);
    assert!(output.stdout == known_out.as_bytes().to_vec());
}
