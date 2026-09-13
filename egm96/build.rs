use std::fs::{self, File};
use std::io::Write;
use std::path::Path;

const URL_ROOT: &str = "https://micahcc.github.io/egm96-rs/egm96/data";

#[allow(unused)]
fn load_blob(name: &str, env_name: &str, url: String, out_name: String) {
    if let Ok(env) = std::env::var(env_name) {
        std::fs::copy(&env, &out_name)
            .unwrap_or_else(|_| panic!("Failed to copy file from {env} to {out_name}"));
        return;
    }

    #[cfg(feature = "fetch-maps")]
    {
        let help = format!("To use a local file: set environment variable: {env_name}");

        let response = reqwest::blocking::get(&url)
            .unwrap_or_else(|_| panic!("Failed to GET {name}. {help}"))
            .error_for_status()
            .unwrap_or_else(|_| panic!("Failed HTTP status for {name}. {help}"));

        let content = response
            .bytes()
            .unwrap_or_else(|_| panic!("Failed to read bytes for {name}. {help}"));

        let mut dest = File::create(&out_name)
            .unwrap_or_else(|_| panic!("Failed to create output file {out_name}. {help}"));

        dest.write_all(&content)
            .unwrap_or_else(|_| panic!("Failed to write {out_name}. {help}"));

        return;
    }

    panic!("fetch-maps feature is not enabled, and environment {env_name} has not been set!");
}

fn main() {
    println!("cargo:rerun-if-changed=data/coefficients.txt");

    let out_dir = std::env::var("OUT_DIR").unwrap_or_else(|_| panic!("OUT_DIR not set"));

    // 1. Generate coefficients file for egm96_data.rs
    let dest_path = Path::new(&out_dir).join("generated_coefficients.rs");
    let content = fs::read_to_string("data/coefficients.txt")
        .expect("Failed to read coefficients source file");

    let generated = format!(
        "pub static EGM96_DATA: [[f64; 4]; 65342] = [\n{}\n];",
        content
    );

    fs::write(&dest_path, generated).expect("Failed to write generated coefficients");

    // 2. Fetch or copy raster fixtures only for enabled raster features.
    #[cfg(feature = "raster_15_min")]
    load_blob(
        "egm96-15.png",
        "EGM96_15_MIN",
        format!("{URL_ROOT}/egm96-15.png"),
        format!("{out_dir}/egm96-15.png"),
    );

    #[cfg(feature = "raster_5_min")]
    load_blob(
        "egm96-5.png",
        "EGM96_5_MIN",
        format!("{URL_ROOT}/egm96-5.png"),
        format!("{out_dir}/egm96-5.png"),
    );
}
