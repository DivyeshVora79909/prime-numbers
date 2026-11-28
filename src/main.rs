use clap::{Parser, Subcommand};
use prime_label::{Changes, Integer, PrimeId};
use std::io::{self, BufRead};
use std::path::PathBuf;
use std::process::ExitCode;

#[derive(Parser)]
#[command(version, about)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Generate a random, proven 40-digit prime ID.
    NewId,
    /// Create the empty set without replacing an existing file.
    Init { label: PathBuf },
    /// Print true or false for a prime ID in this label.
    Contains { label: PathBuf, id: PrimeId },
    /// Print the decimal remainder for a positive modulus.
    Mod { label: PathBuf, modulus: Integer },
    /// Apply +prime grants and -prime revocations together. Without arguments, read changes from stdin.
    Apply {
        label: PathBuf,
        #[arg(value_name = "CHANGE", allow_hyphen_values = true)]
        changes: Vec<String>,
    },
}

fn request(batch: &mut Changes, token: &str) -> io::Result<()> {
    if let Some(id) = token.strip_prefix('+') {
        batch.grant(id.parse()?)
    } else if let Some(id) = token.strip_prefix('-') {
        batch.revoke(id.parse()?)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "expected +prime or -prime",
        ))
    }
}

fn run(command: Command) -> io::Result<()> {
    match command {
        Command::NewId => println!("{}", PrimeId::generate()?),
        Command::Init { label } => prime_label::init(label)?,
        Command::Contains { label, id } => println!("{}", prime_label::contains(label, &id)?),
        Command::Mod { label, modulus } => println!("{}", prime_label::remainder(label, &modulus)?),
        Command::Apply { label, changes } => {
            let mut batch = Changes::default();
            if changes.is_empty() {
                for line in io::stdin().lock().lines() {
                    for token in line?.split_ascii_whitespace() {
                        request(&mut batch, token)?;
                    }
                }
            } else {
                for token in changes {
                    request(&mut batch, &token)?;
                }
            }
            let changed = prime_label::apply(label, batch)?;
            println!("{}", if changed { "changed" } else { "unchanged" });
        }
    }
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse().command) {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("prime-label: {err}");
            ExitCode::FAILURE
        }
    }
}
