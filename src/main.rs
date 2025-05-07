use clap::Parser;
use clap_verbosity_flag::{InfoLevel, Verbosity};
use inversion_finder::*;
use log::info;
use std::path::PathBuf;

/// Look for inversions in a pangenome graph in GFA format
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// path to input gfa
    gfa: PathBuf,

    /// name of reference path
    ref_path: String,

    /// when aligning paths longer than this, use the lowmem algorithm
    #[arg(short, long, default_value_t = 10000)]
    max_highmem_path_length: usize,

    /// maximum path length to align
    #[arg(short = 'p', long, default_value_t = 100000)]
    max_path_length: usize,

    /// minimum length of an inversion in bp for it to be reported
    #[arg(short = 'l', long, default_value_t = 50)]
    min_inversion_length: i32,

    /// maximum drop for heuristic in lowmem mode
    #[arg(short = 'd', long, default_value_t = 1000)]
    max_lowmem_drop: usize,

    /// comma-separated list of paths to exclude
    #[arg(short, long, default_value = "")]
    exclude: String,

    /// verbosity level of logging to stderr
    #[command(flatten)]
    verbose: Verbosity<InfoLevel>,
}

fn main() {
    let args = Args::parse();

    stderrlog::new()
        .module(module_path!())
        .verbosity(args.verbose.log_level_filter())
        .timestamp(stderrlog::Timestamp::Second)
        .init()
        .unwrap();

    info!("Reading GFA");
    let (segment_lengths, paths, path_names) = gfa::read_gfa(args.gfa);

    let ref_path_key = if paths.contains_key(&args.ref_path) {
        args.ref_path
    } else {
        paths
            .keys()
            .find(|k| k.split("#").collect::<Vec<_>>()[0] == args.ref_path)
            .unwrap()
            .clone()
    };

    let paths_to_exclude: Vec<_> = args.exclude.split(",").collect();

    let (inversions, query_path_keys) = alignment_interface::align_all_queries(
        &segment_lengths,
        &paths,
        &path_names,
        &paths_to_exclude,
        &ref_path_key,
        alignment_interface::AlignmentOptions {
            max_highmem_path_length: args.max_highmem_path_length,
            max_lowmem_drop: args.max_lowmem_drop,
            max_path_length: args.max_path_length,
        },
    );

    alignment_interface::print_collated_inversions(
        &inversions,
        &query_path_keys,
        &ref_path_key,
        args.min_inversion_length,
    );
}
