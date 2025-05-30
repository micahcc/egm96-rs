use std::fs::File;
use std::io::Write;

const URL_ROOT: &str = "https://micahcc.github.io/egm96-rs/egm96/data";

struct Fixture<'a> {
    name: &'a str,
    environ: &'a str,
}

fn load_blob(name: &str, env_name: &str, url: String, out_name: String) {
    if let Ok(env) = std::env::var(env_name) {
        std::os::unix::fs::symlink(env, out_name).expect("Failed to symlink.");
    } else {
        let help = format!("To use a local file: set environment variable: {env_name}");
        let response =
            reqwest::blocking::get(url).expect(&format!("Failed to get {}. {help}", name));
        let content = response.bytes().expect("Failed to get bytes. {help}");
        let mut dest =
            File::create(&out_name).expect(&format!("Failed to create output file. {help}"));
        dest.write_all(&content)
            .expect(&format!("Failed to write {out_name}. {help}"));
    }
}

fn main() {
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
            format!(
                "{}/{}",
                std::env::var("OUT_DIR").expect("no OUT_DIR"),
                fixture.name
            ),
        );
    }
}
