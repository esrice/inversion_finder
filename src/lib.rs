use ndarray::{Array, Array2};
use ndarray_stats::QuantileExt;
use regex::Regex;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

pub fn lookup_base_positions(
    path: &[i32],
    segment_lengths: &HashMap<i32, i32>,
    segments: &[i32],
) -> HashMap<i32, (i32, i32)> {
    // first, populate the map of segment positions with the segments we actually want to lookup as
    // keys and dummy values
    let mut segment_positions: HashMap<i32, (i32, i32)> = HashMap::new();
    for segment in segments {
        segment_positions.insert(*segment, (-1, -1));
    }

    let mut current_position = 0;
    for segment in path {
        let this_segment_length = segment_lengths.get(segment).unwrap();
        if segment_positions.contains_key(segment) {
            segment_positions.insert(
                *segment,
                (current_position + 1, current_position + this_segment_length),
            );
        }
        current_position += this_segment_length;
    }
    return segment_positions;
}

/// Parse the path part of a GFA P-line.
///
/// # Arguments
///
/// * `path_string` - the path field from a P-line of a GFA
///
/// # Returns
///
/// The path as a vector of segment IDs, where the orientation is encoded by the sign of the
/// segment ID. So, for example, "1-" is encoded as -1, while "1+" is encoded as 1.
///
/// # Examples
///
/// ```
/// assert_eq!(inversion_finder::parse_gfa_path("1-,2+,3-,235+"), vec![-1, 2, -3, 235]);
/// ```
pub fn parse_gfa_path(path_string: &str) -> Vec<i32> {
    let re = Regex::new(r"(\d+)([+-])").unwrap();
    let mut path_list: Vec<i32> = Vec::new();
    for segment in path_string.split(",") {
        let caps = re.captures(segment).unwrap();
        path_list.push(match caps.get(2).unwrap().as_str() {
            "+" => caps.get(1).unwrap().as_str().parse::<i32>().unwrap(),
            "-" => -1 * caps.get(1).unwrap().as_str().parse::<i32>().unwrap(),
            _ => panic!("Invalid segment {}", segment),
        })
    }
    return path_list;
}

/// Read a GFA into memory.
///
/// Only keep the information in the GFA that we will need later: the length of each segment, and
/// the paths.
pub fn read_gfa(path: PathBuf) -> (HashMap<i32, i32>, HashMap<String, Vec<i32>>) {
    let file = match File::open(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.display(), why),
        Ok(file) => file,
    };
    let reader = BufReader::new(file);

    let mut segment_lengths: HashMap<i32, i32> = HashMap::new();
    let mut paths = HashMap::new();

    for line in reader.lines().map(|l| l.unwrap()) {
        let fields: Vec<&str> = line.split("\t").collect();

        if fields[0] == "S" {
            segment_lengths.insert(
                fields[1].to_string().parse().unwrap(),
                fields[2].len().try_into().unwrap(),
            );
        } else if fields[0] == "P" {
            paths.insert(fields[1].to_string(), parse_gfa_path(fields[2]));
        }
    }

    return (segment_lengths, paths);
}

// This is a janky function that should only be used for finding the max of an array of possible
// scores and for no other purpose.
fn score_max(a: [i32; 4]) -> i32 {
    *a.iter().max().unwrap()
}

// This is a janky function that should only be used for finding the argmax of an array of possible
// scores and for no other purpose.
fn score_argmax(a: [i32; 4]) -> i8 {
    (0..4).max_by_key(|x| a[*x]).unwrap().try_into().unwrap()
}

/// Create alignment matrices with edges filled.
pub fn create_matrices(
    path1: &[i32],
    path2: &[i32],
    segment_lengths: &HashMap<i32, i32>,
) -> (Array2<i32>, Array2<i8>) {
    let mut score_matrix: Array2<i32> = Array::zeros((path1.len(), path2.len()));
    let mut traceback_matrix: Array2<i8> = Array::zeros((path1.len(), path2.len()));

    // fill in the corner
    score_matrix[[0, 0]] = if path1[0] == path2[0] {
        *segment_lengths.get(&path1[0].abs()).unwrap()
    } else {
        -1
    };

    // fill in the first column
    for i in 1..path1.len() {
        let this_cell_score = if path1[i] == path2[0] {
            *segment_lengths.get(&path1[i]).unwrap()
        } else {
            -1
        };

        let possible_scores = [0, -1, score_matrix[[i - 1, 0]], -1];
        score_matrix[[i, 0]] = score_max(possible_scores) + this_cell_score;
        traceback_matrix[[i, 0]] = score_argmax(possible_scores);
    }

    // fill in the first row
    for j in 1..path2.len() {
        let this_cell_score = if path2[j] == path1[0] {
            *segment_lengths.get(&path2[j]).unwrap()
        } else {
            -1
        };

        let possible_scores = [0, -1, -1, score_matrix[[0, j - 1]]];
        score_matrix[[0, j]] = score_max(possible_scores) + this_cell_score;
        traceback_matrix[[0, j]] = score_argmax(possible_scores);
    }

    return (score_matrix, traceback_matrix);
}

/// Perform an alignment subproblem.
///
/// An alignment subproblem is one where neither of the paths contains a segment that is present in
/// opposite orientations in the two paths.
fn align_paths_subproblem(
    path1: &[i32],
    path2: &[i32],
    segment_lengths: &HashMap<i32, i32>,
) -> (Vec<i32>, Vec<i32>) {
    let (mut score_matrix, mut traceback_matrix) =
        create_matrices(&path1, &path2, &segment_lengths);

    for i in 1..path1.len() {
        for j in 1..path2.len() {
            let this_cell_score = if path1[i] == path2[j] {
                *segment_lengths.get(&path1[i].abs()).unwrap()
            } else {
                -1
            };

            let possible_scores = [
                0,
                score_matrix[[i - 1, j - 1]], // come from diagonal
                score_matrix[[i - 1, j]],     // come from above
                score_matrix[[i, j - 1]],     // come from left
            ];

            score_matrix[[i, j]] = score_max(possible_scores) + this_cell_score;
            traceback_matrix[[i, j]] = score_argmax(possible_scores);
        }
    }
    return traceback(&path1, &path2, &score_matrix, &traceback_matrix);
}

fn traceback(
    path1: &[i32],
    path2: &[i32],
    score_matrix: &Array2<i32>,
    traceback_matrix: &Array2<i8>,
) -> (Vec<i32>, Vec<i32>) {
    let (mut i, mut j) = score_matrix.argmax().unwrap();
    let mut alignment_end_reached = false;
    let mut alignment_path1: Vec<i32> = Vec::new();
    let mut alignment_path2: Vec<i32> = Vec::new();
    while !alignment_end_reached {
        let segment_path1 = path1[i];
        let segment_path2 = path2[j];
        if alignment_path1.is_empty() || *alignment_path1.last().unwrap() != segment_path1 {
            alignment_path1.push(segment_path1);
        }
        if alignment_path2.is_empty() || *alignment_path2.last().unwrap() != segment_path2 {
            alignment_path2.push(segment_path2);
        }

        match traceback_matrix[[i, j]] {
            0 => alignment_end_reached = true,
            1 => {
                i -= 1;
                j -= 1;
            }
            2 => i -= 1,
            3 => j -= 1,
            x => panic!("Bad value {} in traceback matrix!", x),
        }
    }
    alignment_path1.reverse();
    alignment_path2.reverse();
    return (alignment_path1, alignment_path2);
}

pub fn align_paths(
    path1: &[i32],
    path2: &[i32],
    segment_lengths: &HashMap<i32, i32>,
) -> Vec<(Vec<i32>, Vec<i32>)> {
    // reverse-complemented version of path2
    let path2_rev: Vec<i32> = path2.iter().map(|x| -1 * x).rev().collect();
    let path1_set = HashSet::<_>::from_iter(path1.iter().cloned());
    let path2_set = HashSet::<_>::from_iter(path2.iter().cloned());
    let path2_rev_set = HashSet::<_>::from_iter(path2_rev.iter().cloned());
    // all segments traversed in the same direction by path1 and path2
    let conflicting_segments =
        HashSet::<_>::from_iter(path1_set.intersection(&path2_set).map(|x| x.abs()));
    // all segments traversed in opposite directions by path1 and path2, and NOT also traversed in
    // the same direction by the two paths somewhere else
    let intersection =
        HashSet::<_>::from_iter(path1_set.intersection(&path2_rev_set).map(|x| x.abs()));
    let common_segments = HashSet::<_>::from_iter(intersection.difference(&conflicting_segments));
    let mut used_segments: HashSet<i32> = HashSet::new();

    let mut alignments = Vec::new();

    for i in 0..path1.len() {
        if common_segments.contains(&path1[i].abs()) && !used_segments.contains(&path1[i].abs()) {
            let mut k = i;
            while k < path1.len() && !conflicting_segments.contains(&path1[k].abs()) {
                k += 1;
            }
            let path1_subproblem = &path1[i..k];

            let j = path2_rev.iter().position(|&x| x == path1[i]).unwrap();
            k = j;
            while k < path2_rev.len() && !conflicting_segments.contains(&path2_rev[k].abs()) {
                k += 1;
            }
            let path2_subproblem = &path2_rev[j..k];

            let (alignment_path1, alignment_path2) =
                align_paths_subproblem(path1_subproblem, path2_subproblem, segment_lengths);
            for segment in &alignment_path1 {
                used_segments.insert(segment.abs());
            }
            for segment in &alignment_path2 {
                used_segments.insert(segment.abs());
            }
            alignments.push((
                alignment_path1,
                alignment_path2.iter().rev().map(|x| -1 * x).collect(),
            ));
        }
    }

    return alignments;
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    #[should_panic]
    fn test_parse_bad_gfa_path() {
        parse_gfa_path("1-,2+,3-,235");
    }

    #[test]
    fn test_create_matrices() {
        let path1 = vec![2, 3, 4, -5];
        let path2 = vec![2, 7, -5];
        let segment_lengths: HashMap<i32, i32> =
            HashMap::from([(2, 100), (3, 100), (4, 100), (5, 100), (7, 100)]);
        let (score_matrix, traceback_matrix) = create_matrices(&path1, &path2, &segment_lengths);
        assert_eq!(
            score_matrix,
            array![[100, 99, 98], [99, 0, 0], [98, 0, 0], [97, 0, 0]]
        );

        assert_eq!(
            traceback_matrix,
            array![[0, 3, 3], [2, 0, 0], [2, 0, 0], [2, 0, 0]],
        );
    }

    #[test]
    fn test_align_subproblem() {
        let path1 = vec![2, 3, 4, -5, 6];
        let path2 = vec![6, 2, 7, -5];
        let mut segment_lengths: HashMap<i32, i32> = HashMap::new();
        for i in 0..7 {
            segment_lengths.insert(i, 100);
        }
        let (path1_alignment, path2_alignment) =
            align_paths_subproblem(&path1, &path2, &segment_lengths);
        assert_eq!(path1_alignment, vec![2, 3, 4, -5]);
        assert_eq!(path2_alignment, vec![2, 7, -5]);
    }

    #[test]
    fn test_align_paths() {
        let path1 = vec![1, 2, 3, 4, 5, 6];
        let path2 = vec![1, -5, -7, -2, 6];
        let mut segment_lengths: HashMap<i32, i32> = HashMap::new();
        for i in 0..10 {
            segment_lengths.insert(i, 100);
        }

        assert_eq!(
            align_paths(&path1, &path2, &segment_lengths),
            vec![(vec![2, 3, 4, 5], vec![-5, -7, -2])]
        );

        let path3 = vec![1, 2, 3, 4, 5, 6, 7];
        let path4 = vec![1, -3, -2, 4, -6, -5, 7];
        assert_eq!(
            align_paths(&path3, &path4, &segment_lengths),
            vec![(vec![2, 3], vec![-3, -2]), (vec![5, 6], vec![-6, -5])]
        );

        let path5 = vec![1, 2, 3, 4, 5, 6, 7];
        let path6 = vec![1, -3, -2, 8, -6, -5, 7];
        assert_eq!(
            align_paths(&path5, &path6, &segment_lengths),
            vec![(vec![2, 3], vec![-3, -2]), (vec![5, 6], vec![-6, -5])]
        );
    }
}
