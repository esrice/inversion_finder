extern crate clap;

use std::collections::HashMap;
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

    #[arg(short, long)]
    ref_path: String,
}

fn main() {
    let args = Args::parse();

    let (segment_lengths, paths) = read_gfa(args.gfa);

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

    let mut query_path_keys = Vec::<String>::new();
    let mut inversions = Vec::<(String, i32, i32)>::new();
    for (query_path_key, query_path) in paths {
        if query_path_key != ref_path_key {
            query_path_keys.push(query_path_key.clone());
            let alignments = align_paths(&ref_path, &query_path, &segment_lengths);

            // make a list of segments that we need to find the positions of
            let mut segments_to_lookup = Vec::new();
            for alignment in &alignments {
                segments_to_lookup.push(alignment.0[0]);
                segments_to_lookup.push(alignment.0[alignment.0.len() - 1]);
            }
            let base_positions =
                lookup_base_positions(&ref_path, &segment_lengths, &segments_to_lookup);

            for alignment in alignments {
                inversions.push((
                    query_path_key.clone(),
                    base_positions.get(&alignment.0[0]).unwrap().0,
                    base_positions
                        .get(&alignment.0[alignment.0.len() - 1])
                        .unwrap()
                        .1,
                ));
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

    println!("ref\tstart\tend\t{}", query_path_keys.join("\t"));
    for ((start_position, end_position), paths) in inversions_collated {
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
