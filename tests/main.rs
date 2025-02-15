use assert_cmd::prelude::*;
use std::process::Command;

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
        .args(&["-i", "example_data/test.mtx", "--atac", "example_data/target.bed.gz", "--barcodes", "example_data/test.tsv", "--nchr", "1", "--loglevel", "warn", "--bed", "example_data/test.bed", "--cores", "1"])
        .assert();

    let output = out
        .get_output();
    let known_out = "barcode1 0.33333\nbarcode2 0.14286\n";
    let known_out2 = "barcode2 0.14286\nbarcode1 0.33333\n";

    // Seems to be somewhat random which order
    // the cells are printed in, so just check for
    // either.
    assert!(output.stdout == known_out.as_bytes().to_vec() || output.stdout == known_out2.as_bytes().to_vec() );
}

#[test]
fn test_full_output() {
 let out = sc_jaccard()
        .args(&["-i", "example_data/test.mtx", "--atac", "example_data/target.bed.gz", "--barcodes", "example_data/test.tsv", "--nchr", "1", "--loglevel", "warn", "--bed", "example_data/test.bed", "--cores", "1", "--full"])
        .assert();

    let output = out
        .get_output();
    let known_out = "barcode1, 200, 100, 2, 200, 2, 0.33333\nbarcode2, 200, 50, 1, 200, 2, 0.14286\n";
    let known_out2 = "barcode2, 200, 50, 1, 200, 2, 0.14286\nbarcode1, 200, 100, 2, 200, 2, 0.33333\n";

    // Seems to be somewhat random which order
    // the cells are printed in, so just check for
    // either.
    assert!(output.stdout == known_out.as_bytes().to_vec() || output.stdout == known_out2.as_bytes().to_vec() );

}
