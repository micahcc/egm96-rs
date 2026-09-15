use assert_cmd::Command;
use assert_fs::prelude::*;
use predicates::prelude::*;

#[test]
fn test_cli_single_point_flags() {
    let mut cmd = Command::cargo_bin("egm96-cli").unwrap();
    cmd.arg("--lat")
        .arg("37.7749")
        .arg("--lon")
        .arg("-122.4194")
        .assert()
        .success()
        .stdout(predicate::str::contains("-32.1904"));
}

#[test]
fn test_cli_single_point_subcommand() {
    let mut cmd = Command::cargo_bin("egm96-cli").unwrap();
    cmd.arg("point")
        .arg("37.7749")
        .arg("-122.4194")
        .assert()
        .success()
        .stdout(predicate::str::contains("-32.1904"));
}

#[test]
fn test_cli_missing_args_fails() {
    let mut cmd = Command::cargo_bin("egm96-cli").unwrap();
    cmd.arg("--lat")
        .arg("37.7749")
        .assert()
        .failure()
        .stderr(predicate::str::contains("Error"));
}

#[test]
fn test_cli_lon_only_fails() {
    let mut cmd = Command::cargo_bin("egm96-cli").unwrap();
    cmd.arg("--lon")
        .arg("-122.4194")
        .assert()
        .failure()
        .stderr(predicate::str::contains("Error"));
}

#[test]
fn test_cli_batch_stdin() {
    let mut cmd = Command::cargo_bin("egm96-cli").unwrap();
    cmd.arg("batch")
        .write_stdin("latitude,longitude\n37.7749,-122.4194\n40.7128,-74.0060\n")
        .assert()
        .success()
        .stdout(predicate::str::contains("37.7749,-122.4194,-32.1904"))
        .stdout(predicate::str::contains("40.7128,-74.006"));
}

#[test]
fn test_cli_batch_file_io() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let input_file = temp_dir.child("coords.csv");
    input_file
        .write_str("latitude,longitude\n37.7749,-122.4194\n")
        .unwrap();

    let output_file = temp_dir.child("results.csv");

    let mut cmd = Command::cargo_bin("egm96-cli").unwrap();
    cmd.arg("batch")
        .arg("--input")
        .arg(input_file.path())
        .arg("--output")
        .arg(output_file.path())
        .assert()
        .success();

    output_file.assert(predicate::path::exists());
    output_file.assert(predicate::str::contains("37.7749,-122.4194,-32.1904"));
}

#[test]
fn test_cli_batch_input_only() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let input_file = temp_dir.child("coords.csv");
    input_file
        .write_str("latitude,longitude\n37.7749,-122.4194\n")
        .unwrap();

    let mut cmd = Command::cargo_bin("egm96-cli").unwrap();
    cmd.arg("batch")
        .arg("--input")
        .arg(input_file.path())
        .assert()
        .success()
        .stdout(predicate::str::contains("37.7749,-122.4194,-32.1904"));
}

#[test]
fn test_cli_batch_output_only() {
    let temp_dir = assert_fs::TempDir::new().unwrap();
    let output_file = temp_dir.child("results.csv");

    let mut cmd = Command::cargo_bin("egm96-cli").unwrap();
    cmd.arg("batch")
        .arg("--output")
        .arg(output_file.path())
        .write_stdin("latitude,longitude\n37.7749,-122.4194\n")
        .assert()
        .success();

    output_file.assert(predicate::path::exists());
    output_file.assert(predicate::str::contains("37.7749,-122.4194,-32.1904"));
}

#[test]
fn test_cli_batch_nonexistent_input() {
    let mut cmd = Command::cargo_bin("egm96-cli").unwrap();
    cmd.arg("batch")
        .arg("--input")
        .arg("nonexistent_file_path_12345.csv")
        .assert()
        .failure();
}

#[test]
fn test_cli_batch_malformed_csv_skipped() {
    let mut cmd = Command::cargo_bin("egm96-cli").unwrap();
    cmd.arg("batch")
        .write_stdin("latitude,longitude\ninvalid_row\n37.7749,-122.4194\n")
        .assert()
        .success()
        .stdout(predicate::str::contains("37.7749,-122.4194,-32.1904"));
}

#[test]
fn test_cli_batch_extra_columns_handled() {
    let mut cmd = Command::cargo_bin("egm96-cli").unwrap();
    cmd.arg("batch")
        .write_stdin("latitude,longitude\n37.7749,-122.4194,extra_data\n40.7128,-74.0060\n")
        .assert()
        .success()
        .stdout(predicate::str::contains("37.7749,-122.4194,-32.1904"))
        .stdout(predicate::str::contains("40.7128,-74.006"));
}

#[test]
fn test_cli_batch_parse_error_fails() {
    let mut cmd = Command::cargo_bin("egm96-cli").unwrap();
    cmd.arg("batch")
        .write_stdin("latitude,longitude\nabc,xyz\n")
        .assert()
        .failure();
}
