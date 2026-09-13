use std::env;
use std::path::Path;

fn main() {
    // just compile with the EGM96.h to ensure it is valid c
    println!("cargo:rerun-if-changed=src/lib.rs");
    println!("cargo:rerun-if-changed=Cargo.toml");

    let manifest_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
    let header_path = Path::new(&manifest_dir).join("src/EGM96.h");

    let config = cbindgen::Config {
        pragma_once: true,
        language: cbindgen::Language::C,
        cpp_compat: true,
        tab_width: 2,
        braces: cbindgen::Braces::NextLine,
        ..Default::default()
    };

    let result = cbindgen::Builder::new()
        .with_crate(&manifest_dir)
        .with_config(config)
        .with_pragma_once(true)
        .generate();

    match result {
        Ok(binding) => {
            binding.write_to_file(header_path);
        }
        Err(err) => {
            panic!("Failed to generate C bindings via cbindgen: {err}");
        }
    }
}
