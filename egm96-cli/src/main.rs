use clap::{Parser, Subcommand};
use std::fs::File;
use std::io;
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    name = "egm96-cli",
    about = "Query EGM96 geoid undulation height offsets."
)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,

    /// Latitude shortcut
    #[arg(short, long, allow_hyphen_values = true)]
    lat: Option<f64>,

    /// Longitude shortcut
    #[arg(short = 'o', long, allow_hyphen_values = true)]
    lon: Option<f64>,
}

#[derive(Subcommand)]
enum Commands {
    Point {
        #[arg(allow_hyphen_values = true)]
        lat: f64,
        #[arg(allow_hyphen_values = true)]
        lon: f64,
    },
    Batch {
        #[arg(short, long)]
        input: Option<PathBuf>,
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    match (cli.command, cli.lat, cli.lon) {
        (None, Some(lat), Some(lon)) | (Some(Commands::Point { lat, lon }), _, _) => {
            let offset = egm96::egm96_altitude_offset(lat, lon);
            println!("{:.4}", offset);
        }
        (Some(Commands::Batch { input, output }), _, _) => {
            process_batch(input, output)?;
        }
        _ => {
            eprintln!("Error: Provide --lat and --lon, or use 'point'/'batch' subcommands.");
            std::process::exit(1);
        }
    }
    Ok(())
}

fn process_batch(
    input_path: Option<PathBuf>,
    output_path: Option<PathBuf>,
) -> Result<(), Box<dyn std::error::Error>> {
    let input_reader: Box<dyn io::Read> = match input_path {
        Some(path) => Box::new(File::open(path)?),
        None => Box::new(io::stdin()),
    };
    let mut rdr = csv::ReaderBuilder::new()
        .flexible(true)
        .from_reader(input_reader);

    let output_writer: Box<dyn io::Write> = match output_path {
        Some(path) => Box::new(File::create(path)?),
        None => Box::new(io::stdout()),
    };
    let mut wtr = csv::Writer::from_writer(output_writer);

    wtr.write_record(["latitude", "longitude", "offset_m"])?;

    for result in rdr.records() {
        let record = result?;
        if record.len() >= 2 {
            let lat: f64 = record[0].parse()?;
            let lon: f64 = record[1].parse()?;
            let offset = egm96::egm96_altitude_offset(lat, lon);
            wtr.write_record([lat.to_string(), lon.to_string(), format!("{:.4}", offset)])?;
        }
    }

    wtr.flush()?;
    Ok(())
}
