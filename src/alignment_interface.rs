use std::{collections::HashMap, error::Error};

use log::info;

use crate::{InversionError, align, gfa};

pub struct AlignmentOptions {
    pub max_highmem_path_length: usize,
    pub max_lowmem_drop: usize,
    pub max_path_length: usize,
}

/// Align every non-reference path to the reference.
///
/// # Arguments
///
/// * `segment_lengths`: a map of segment ID to segment length in bp
/// * `paths`: map of path name to path, represented as sequence of path IDs
/// * `path_names`: keys of `paths`, in the order they should be accessed
/// * `paths_to_exclude`: keys of all paths that should not be aligned
/// * `ref_path_key`: key in `paths` of the reference path
/// * `alignment_options`: parameters for the alignments
///
/// # Returns
///
/// * `inversions`: a vec of tuples (query path key, node ID of first node in inverison, node ID of
///   last node in inversion)
/// * `query_path_keys`: a vec of keys for paths which were actually aligned to the reference
pub fn align_all_queries(
    segment_lengths: &HashMap<i32, i32>,
    paths: &HashMap<String, Vec<i32>>,
    path_names: &[String],
    paths_to_exclude: &[&str],
    ref_path_key: &str,
    alignment_options: AlignmentOptions,
) -> Result<(Vec<(String, i32, i32)>, Vec<String>), InversionError> {
    let mut query_path_keys = Vec::<String>::new();
    let mut inversions = Vec::<(String, i32, i32)>::new();
    let ref_path = paths
        .get(ref_path_key)
        .ok_or(InversionError::PathNotFound(ref_path_key.to_string()))?
        .clone();

    for query_path_key in path_names {
        if query_path_key != ref_path_key
            && !paths_to_exclude
                .iter()
                .any(|x| x == query_path_key || *x == query_path_key.split("#").nth(0).unwrap())
        {
            info!("Starting alignment of path {}", query_path_key);
            let query_path = paths
                .get(query_path_key)
                .ok_or(InversionError::PathNotFound(query_path_key.to_string()))?;
            query_path_keys.push(query_path_key.clone());
            let alignments = align::align_paths(
                &ref_path,
                &query_path,
                &segment_lengths,
                alignment_options.max_highmem_path_length,
                alignment_options.max_lowmem_drop,
                alignment_options.max_path_length,
            )?;

            // make a list of segments that we need to find the positions of
            let mut segments_to_lookup = Vec::new();
            for alignment in &alignments {
                segments_to_lookup.push(alignment.path1_start_index);
                segments_to_lookup.push(alignment.path1_end_index);
            }
            let base_positions =
                gfa::lookup_base_positions(&ref_path, &segment_lengths, &segments_to_lookup)?;

            for alignment in alignments {
                let start_position = base_positions.get(&alignment.path1_start_index).unwrap().0;
                let end_position = base_positions.get(&alignment.path1_end_index).unwrap().1;
                inversions.push((query_path_key.clone(), start_position, end_position));
            }
        }
    }
    Ok((inversions, query_path_keys))
}

pub fn print_collated_inversions(
    inversions: &[(String, i32, i32)],
    query_path_keys: &[String],
    ref_path_key: &str,
    min_inversion_length: i32,
) -> Result<(), Box<dyn Error>> {
    // finally, collate the inversions from the different animals and print out a table
    let mut inversions_collated: HashMap<(i32, i32), Vec<String>> = HashMap::new();
    for (path, start_position, end_position) in inversions {
        inversions_collated
            .entry((*start_position, *end_position))
            .and_modify(|v| v.push(path.clone()))
            .or_insert(vec![path.clone()]);
    }

    // print the collated inversions out ordered by start position
    println!("ref\tstart\tend\t{}", query_path_keys.join("\t"));
    let mut keys: Vec<&(i32, i32)> = inversions_collated.keys().collect();
    keys.sort_by_key(|k| k.0);
    for (start_position, end_position) in keys {
        if end_position - start_position >= min_inversion_length {
            let paths = inversions_collated
                .get(&(*start_position, *end_position))
                .ok_or(format!(
                    "Cannot find inversion {}-{}",
                    *start_position, *end_position
                ))?;
            let mut calls = Vec::<i8>::new();
            for query_path_key in query_path_keys {
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
    Ok(())
}
