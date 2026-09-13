use std::fs::{self, File};
use std::io::Write;
use std::path::Path;

const URL_ROOT: &str = "https://micahcc.github.io/egm96-rs/egm96/data";

struct Fixture<'a> {
    name: &'a str,
    environ: &'a str,
}

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

        let response =
            reqwest::blocking::get(&url).unwrap_or_else(|_| panic!("Failed to GET {name}. {help}"));

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

    // 2. Fetch or copy test map fixtures
    let fixtures = [
        Fixture {
            name: "egm96-15.png",
            environ: "EGM96_15_MIN",
        },
        Fixture {
            name: "egm96-5.png",
            environ: "EGM96_5_MIN",
        },
    ];

    for fixture in fixtures {
        load_blob(
            fixture.name,
            fixture.environ,
            format!("{URL_ROOT}/{}", fixture.name),
            format!("{out_dir}/{}", fixture.name),
        );
    }
}
