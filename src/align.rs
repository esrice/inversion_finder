use ndarray::{Array, Array2};
use ndarray_stats::QuantileExt;
use std::collections::{HashMap, HashSet};

use super::{amax, argmax, lowmem};

pub struct Alignment {
    /// segments aligned in path1
    pub alignment_path1: Vec<i32>,

    /// segments aligned in path2
    pub alignment_path2: Vec<i32>,

    /// start index of alignment in path1
    pub path1_start_index: i32,

    /// end index of alignment in path1
    pub path1_end_index: i32,
}

/// Create alignment matrices with edges filled.
///
/// Arguments:
///
/// * `path1` and `path2`: vectors of oriented segments in paths to align
/// * `segment_lengths`: map of segment ID to segment length
///
/// Returns:
///
/// * `score_matrix`: a DP alignment matrix for score with first row and column filled
/// * `traceback_matrix`: a DP alignment matrix for traceback with first row and column filled.
///   Values: 0 => alignment starts here, 1 => alignment comes from diagonal, 2 => alignment comes
///   from above, 3 => alignment comes from left
fn create_matrices(
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
        -1 * (*segment_lengths.get(&path1[0].abs()).unwrap()
            + *segment_lengths.get(&path2[0].abs()).unwrap())
    };

    // fill in the first column
    for i in 1..path1.len() {
        let this_cell_score = if path1[i] == path2[0] {
            *segment_lengths.get(&path1[i].abs()).unwrap()
        } else {
            -1 * *segment_lengths.get(&path1[i].abs()).unwrap()
        };

        let possible_scores = [0, -1, score_matrix[[i - 1, 0]], -1];
        score_matrix[[i, 0]] = amax(&possible_scores) + this_cell_score;
        traceback_matrix[[i, 0]] = argmax(&possible_scores).try_into().unwrap();
    }

    // fill in the first row
    for j in 1..path2.len() {
        let this_cell_score = if path2[j] == path1[0] {
            *segment_lengths.get(&path2[j].abs()).unwrap()
        } else {
            -1 * *segment_lengths.get(&path2[j].abs()).unwrap()
        };

        let possible_scores = [0, -1, -1, score_matrix[[0, j - 1]]];
        score_matrix[[0, j]] = amax(&possible_scores) + this_cell_score;
        traceback_matrix[[0, j]] = argmax(&possible_scores).try_into().unwrap();
    }

    return (score_matrix, traceback_matrix);
}

/// Perform an alignment subproblem.
///
/// An alignment subproblem is one where neither of the paths contains a segment that is present in
/// opposite orientations in the two paths.
///
/// # Arguments
///
/// * `path1` and `path2`: paths to align
/// * `segment_lengths`: map of segment ID to segment length in bp
///
/// # Returns
///
/// * `alignment_path1` and `alignment_path2`: alignment for both paths
/// * `path1_start_index` and `path1_end_index`: indices of start and end segments of alignment in
///    path1
fn align_paths_subproblem(
    path1: &[i32],
    path2: &[i32],
    segment_lengths: &HashMap<i32, i32>,
) -> Alignment {
    let (mut score_matrix, mut traceback_matrix) =
        create_matrices(&path1, &path2, &segment_lengths);

    for i in 1..path1.len() {
        let len_i = *segment_lengths.get(&path1[i].abs()).unwrap();
        for j in 1..path2.len() {
            let len_j = *segment_lengths.get(&path2[j].abs()).unwrap();
            let possible_scores = if path1[i] == path2[j] {
                [
                    len_i,
                    score_matrix[[i - 1, j - 1]] + len_i, // come from diagonal
                    score_matrix[[i - 1, j]] + len_i,     // come from above
                    score_matrix[[i, j - 1]] + len_i,     // come from left
                ]
            } else {
                [
                    -len_i - len_j,
                    score_matrix[[i - 1, j - 1]] - len_i - len_j,
                    score_matrix[[i - 1, j]] - len_i,
                    score_matrix[[i, j - 1]] - len_j,
                ]
            };

            score_matrix[[i, j]] = amax(&possible_scores);
            traceback_matrix[[i, j]] = argmax(&possible_scores).try_into().unwrap();
        }
    }
    return traceback(&path1, &path2, &score_matrix, &traceback_matrix);
}

fn traceback(
    path1: &[i32],
    path2: &[i32],
    score_matrix: &Array2<i32>,
    traceback_matrix: &Array2<i8>,
) -> Alignment {
    let (mut i, mut j) = score_matrix.argmax().unwrap();
    let path1_end_index = i;
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
    let path1_start_index = i;
    alignment_path1.reverse();
    alignment_path2.reverse();
    return Alignment {
        alignment_path1,
        alignment_path2,
        path1_start_index: path1_start_index.try_into().unwrap(),
        path1_end_index: path1_end_index.try_into().unwrap(),
    };
}

pub fn align_paths(
    path1: &[i32],
    path2: &[i32],
    segment_lengths: &HashMap<i32, i32>,
    max_highmem_length: usize,
    max_lowmem_drop: usize,
    max_path_length: usize,
) -> Vec<Alignment> {
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

    for subproblem1_start in 0..path1.len() {
        if common_segments.contains(&path1[subproblem1_start].abs())
            && !used_segments.contains(&path1[subproblem1_start].abs())
        {
            let mut subproblem1_end = subproblem1_start;
            while subproblem1_end < path1.len()
                && !conflicting_segments.contains(&path1[subproblem1_end].abs())
                && !used_segments.contains(&path1[subproblem1_end].abs())
            {
                subproblem1_end += 1;
            }
            let path1_subproblem = &path1[subproblem1_start..subproblem1_end];

            let subproblem2_start = path2_rev
                .iter()
                .position(|&x| x == path1[subproblem1_start])
                .unwrap();
            let mut subproblem2_end = subproblem2_start;
            while subproblem2_end < path2_rev.len()
                && !conflicting_segments.contains(&path2_rev[subproblem2_end].abs())
                && !used_segments.contains(&path2_rev[subproblem2_end].abs())
            {
                subproblem2_end += 1;
            }
            let path2_subproblem = &path2_rev[subproblem2_start..subproblem2_end];

            // choose correct alignment algorithm depending on length
            let alignment_option = if path1_subproblem.len() < max_highmem_length
                && path2_subproblem.len() < max_highmem_length
            {
                Some(align_paths_subproblem(
                    path1_subproblem,
                    path2_subproblem,
                    segment_lengths,
                ))
            } else if path1_subproblem.len() < max_path_length
                && path2_subproblem.len() < max_path_length
            {
                Some(lowmem::align_paths_subproblem_lowmem(
                    path1_subproblem,
                    path2_subproblem,
                    segment_lengths,
                    max_lowmem_drop,
                ))
            } else {
                None
            };

            if let Some(alignment) = alignment_option {
                for segment in &alignment.alignment_path1 {
                    used_segments.insert(segment.abs());
                }
                for segment in &alignment.alignment_path2 {
                    used_segments.insert(segment.abs());
                }
                alignments.push(Alignment {
                    alignment_path1: alignment.alignment_path1,
                    alignment_path2: alignment
                        .alignment_path2
                        .iter()
                        .rev()
                        .map(|x| -1 * x)
                        .collect(),
                    path1_start_index: (subproblem1_start as i32) + alignment.path1_start_index,
                    path1_end_index: (subproblem1_start as i32) + alignment.path1_end_index,
                });
            }
        }
    }

    return alignments;
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn test_create_matrices() {
        let path1 = vec![2, 3, 4, -5];
        let path2 = vec![2, 7, -5];
        let segment_lengths: HashMap<i32, i32> =
            HashMap::from([(2, 100), (3, 10), (4, 10), (5, 100), (7, 10)]);
        let (score_matrix, traceback_matrix) = create_matrices(&path1, &path2, &segment_lengths);
        assert_eq!(
            score_matrix,
            array![[100, 90, -10], [90, 0, 0], [80, 0, 0], [-20, 0, 0]]
        );

        assert_eq!(
            traceback_matrix,
            array![[0, 3, 3], [2, 0, 0], [2, 0, 0], [2, 0, 0]],
        );
    }

    #[test]
    fn test_align_subproblem() {
        let path1: Vec<i32> = vec![2, 3, 4, -5, 6];
        let path2: Vec<i32> = vec![6, 2, 7, -5];
        let mut segment_lengths: HashMap<i32, i32> = HashMap::new();
        for i in 0..8 {
            segment_lengths.insert(
                i,
                if (path1.contains(&i) || path1.contains(&(-i)))
                    && (path2.contains(&i) || path2.contains(&(-i)))
                {
                    100
                } else {
                    10
                },
            );
        }
        //let (path1_alignment, path2_alignment, path1_start_index, path1_end_index) =
        let alignment = align_paths_subproblem(&path1, &path2, &segment_lengths);
        assert_eq!(alignment.alignment_path1, vec![2, 3, 4, -5]);
        assert_eq!(alignment.alignment_path2, vec![2, 7, -5]);
        assert_eq!(alignment.path1_start_index, 0);
        assert_eq!(alignment.path1_end_index, 3);
    }

    #[test]
    fn test_align_paths() {
        let path1 = vec![1, 2, 3, 4, 5, 6];
        let path2 = vec![1, -5, -7, -2, 6];
        let mut segment_lengths: HashMap<i32, i32> = HashMap::new();
        for i in 0..10 {
            segment_lengths.insert(
                i,
                if (path1.contains(&i) || path1.contains(&(-i)))
                    && (path2.contains(&i) || path2.contains(&(-i)))
                {
                    100
                } else {
                    10
                },
            );
        }

        let alignments1 = align_paths(&path1, &path2, &segment_lengths, 10000, 1000, 100000);
        assert_eq!(alignments1[0].alignment_path1, vec![2, 3, 4, 5]);
        assert_eq!(alignments1[0].alignment_path2, vec![-5, -7, -2]);
        assert_eq!(alignments1[0].path1_start_index, 1);
        assert_eq!(alignments1[0].path1_end_index, 4);

        let path3 = vec![1, 2, 3, 4, 5, 6, 7];
        let path4 = vec![1, -3, -2, 4, -6, -5, 7];
        let alignments2 = align_paths(&path3, &path4, &segment_lengths, 10000, 1000, 100000);
        assert_eq!(alignments2[0].alignment_path1, vec![2, 3]);
        assert_eq!(alignments2[0].alignment_path2, vec![-3, -2]);
        assert_eq!(alignments2[0].path1_start_index, 1);
        assert_eq!(alignments2[0].path1_end_index, 2);
        assert_eq!(alignments2[1].alignment_path1, vec![5, 6]);
        assert_eq!(alignments2[1].alignment_path2, vec![-6, -5]);
        assert_eq!(alignments2[1].path1_start_index, 4);
        assert_eq!(alignments2[1].path1_end_index, 5);

        let path5 = vec![1, 2, 3, 4, 5, 6, 7];
        let path6 = vec![1, -3, -2, 8, -6, -5, 7];
        let alignments3 = align_paths(&path5, &path6, &segment_lengths, 10000, 1000, 100000);
        assert_eq!(alignments3[0].alignment_path1, vec![2, 3]);
        assert_eq!(alignments3[0].alignment_path2, vec![-3, -2]);
        assert_eq!(alignments3[0].path1_start_index, 1);
        assert_eq!(alignments3[0].path1_end_index, 2);
        assert_eq!(alignments3[1].alignment_path1, vec![5, 6]);
        assert_eq!(alignments3[1].alignment_path2, vec![-6, -5]);
        assert_eq!(alignments3[1].path1_start_index, 4);
        assert_eq!(alignments3[1].path1_end_index, 5);
    }
}
