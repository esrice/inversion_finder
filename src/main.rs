extern crate clap;

use std::path::PathBuf;

use clap::Parser;

use inversion_finder::*;

/// Look for inversions in a pangenome graph in GFA format
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// path to input gfa
    #[arg(short, long)]
    gfa: PathBuf,
}

fn main() {
    let args = Args::parse();

    read_gfa(args.gfa);
}
