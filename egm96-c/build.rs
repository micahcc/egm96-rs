use std::env;
use std::io::prelude::*;

fn main() {
    // just compile with the EGM96.h to ensure it is valid c
    println!("cargo:rerun-if-changed=src/lib.rs");

    let crate_dir = env::var("CARGO_MANIFEST_DIR").unwrap();

    let config = cbindgen::Config {
        pragma_once: true,
        language: cbindgen::Language::C,
        cpp_compat: true,
        tab_width: 2,
        braces: cbindgen::Braces::NextLine,
        ..Default::default()
    };

    match cbindgen::Builder::new()
        .with_crate(crate_dir)
        .with_config(config)
        .with_pragma_once(true)
        .generate()
    {
        Ok(result) => {
            result.write_to_file("src/EGM96.h");
        }
        Err(err) => {
            let mut file = std::fs::File::create("src/EGM96.h").expect("failed to create");
            file.write_all(format!("/* Failed to compile: {} */", err).as_bytes())
                .expect("failed to write");
        }
    }
}
