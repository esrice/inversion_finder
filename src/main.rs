use clap::Parser;
use clap_verbosity_flag::{InfoLevel, Verbosity};
use inversion_finder::*;
use log::info;
use std::collections::HashMap;
use std::path::PathBuf;

mod gfa;

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

    let ref_path = paths.get(&ref_path_key).unwrap().clone();
    let paths_to_exclude: Vec<_> = args.exclude.split(",").collect();

    let mut query_path_keys = Vec::<String>::new();
    let mut inversions = Vec::<(String, i32, i32)>::new();
    for query_path_key in path_names {
        if query_path_key != ref_path_key
            && !paths_to_exclude
                .iter()
                .any(|x| *x == query_path_key || *x == query_path_key.split("#").nth(0).unwrap())
        {
            info!("Starting alignment of path {}", query_path_key);
            query_path_keys.push(query_path_key.clone());
            let query_path = paths.get(&query_path_key).unwrap();
            let alignments = align::align_paths(
                &ref_path,
                &query_path,
                &segment_lengths,
                args.max_highmem_path_length,
                args.max_lowmem_drop,
            );

            // make a list of segments that we need to find the positions of
            let mut segments_to_lookup = Vec::new();
            for alignment in &alignments {
                segments_to_lookup.push(alignment.path1_start_index);
                segments_to_lookup.push(alignment.path1_end_index);
            }
            let base_positions =
                gfa::lookup_base_positions(&ref_path, &segment_lengths, &segments_to_lookup);

            for alignment in alignments {
                let start_position = base_positions.get(&alignment.path1_start_index).unwrap().0;
                let end_position = base_positions.get(&alignment.path1_end_index).unwrap().1;
                inversions.push((query_path_key.clone(), start_position, end_position));
            }
        }
    }

    // finally, collate the inversions from the different animals and print out a table
    let mut inversions_collated: HashMap<(i32, i32), Vec<String>> = HashMap::new();
    for (path, start_position, end_position) in inversions {
        inversions_collated
            .entry((start_position, end_position))
            .and_modify(|v| v.push(path.clone()))
            .or_insert(vec![path.clone()]);
    }

    // print the collated inversions out ordered by start position
    println!("ref\tstart\tend\t{}", query_path_keys.join("\t"));
    let mut keys: Vec<&(i32, i32)> = inversions_collated.keys().collect();
    keys.sort_by_key(|k| k.0);
    for (start_position, end_position) in keys {
        if end_position - start_position >= args.min_inversion_length {
            let paths = inversions_collated
                .get(&(*start_position, *end_position))
                .unwrap();
            let mut calls = Vec::<i8>::new();
            for query_path_key in &query_path_keys {
                if paths.contains(query_path_key) {
                    calls.push(1);
                } else {
                    calls.push(0);
                }
            }
            println!(
                "{}\t{}\t{}\t{}",
                ref_path_key,
                start_position,
                end_position,
                calls
                    .iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join("\t"),
            );
        }
    }
}
