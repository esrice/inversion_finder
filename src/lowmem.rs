use crate::InversionError;

use super::{align, amax, argmax};
use log::debug;
use std::collections::HashMap;

fn max_and_argmax(a: &[i32]) -> (i32, i32) {
    let mut max = a[0];
    let mut argmax = 0;
    for (i, x) in a.iter().enumerate() {
        if *x > max {
            max = *x;
            argmax = i;
        }
    }
    return (max, argmax.try_into().expect("i32 overflow"));
}

struct InitializeMatricesLowmemResult(Vec<i32>, Vec<i32>, HashMap<(i32, i32), i8>, i32, (i32, i32));

fn initialize_matrices_lowmem(
    path1: &[i32],
    path2: &[i32],
    segment_lengths: &HashMap<i32, i32>,
) -> Result<InitializeMatricesLowmemResult, InversionError> {
    let mut score_row_previous = Vec::<i32>::with_capacity(path2.len());
    // this is a sparse representation of the traceback matrix: if a key isn't present, that
    // indicates that the value is actually 0
    let mut traceback_matrix = HashMap::<(i32, i32), i8>::new();

    // fill in the corner
    score_row_previous.push(if path1[0] == path2[0] {
        *segment_lengths
            .get(&path1[0].abs())
            .ok_or(InversionError::SegmentNotFound(path1[0].abs()))?
    } else {
        -1 * (*segment_lengths
            .get(&path1[0].abs())
            .ok_or(InversionError::SegmentNotFound(path1[0].abs()))?
            + *segment_lengths
                .get(&path2[0].abs())
                .ok_or(InversionError::SegmentNotFound(path2[0].abs()))?)
    });

    // fill in the rest of the first row
    for j in 1..path2.len() {
        let this_cell_score = if path2[j] == path1[0] {
            *segment_lengths
                .get(&path2[j].abs())
                .ok_or(InversionError::SegmentNotFound(path2[j].abs()))?
        } else {
            -1 * *segment_lengths
                .get(&path2[j].abs())
                .ok_or(InversionError::SegmentNotFound(path2[j].abs()))?
        };

        let possible_scores = [0, -1, -1, score_row_previous[j - 1]];
        score_row_previous.push(amax(&possible_scores) + this_cell_score);
        let traceback_value = argmax(&possible_scores).try_into().expect("i32 overflow");
        if traceback_value != 0 {
            traceback_matrix.insert((0, j.try_into().expect("i32 overflow")), traceback_value);
        }
    }

    let (max_score, argmax_score_j) = max_and_argmax(&score_row_previous);
    let argmax_score: (i32, i32) = (0, argmax_score_j);

    let mut score_row_current = Vec::<i32>::with_capacity(path2.len());
    // initialize everything to 0, but these will not actually ever be read. Just for preventing
    // out of bounds errors.
    for _ in 0..path2.len() {
        score_row_current.push(0);
    }

    return Ok(InitializeMatricesLowmemResult(
        score_row_previous,
        score_row_current,
        traceback_matrix,
        max_score,
        argmax_score,
    ));
}

pub fn align_paths_subproblem_lowmem(
    path1: &[i32],
    path2: &[i32],
    segment_lengths: &HashMap<i32, i32>,
    drop: usize,
) -> Result<align::Alignment, InversionError> {
    debug!(
        "Performing lowmem alignment of {}:{} to {}:{} ({}x{})",
        path1[0],
        path1[path1.len() - 1],
        path2[0],
        path2[path2.len() - 1],
        path1.len(),
        path2.len(),
    );
    let InitializeMatricesLowmemResult(
        mut score_row_previous,
        mut score_row_current,
        mut traceback_matrix,
        mut max_score,
        mut argmax_score,
    ) = initialize_matrices_lowmem(&path1, &path2, &segment_lengths)?;

    let (max_row_drop, max_col_drop) = if path1.len() > path2.len() {
        //(min(drop + path1.len() - path2.len(), drop * 5), drop)
        (drop + path1.len() - path2.len(), drop)
    } else {
        //(drop, min(drop + path2.len() - path1.len(), drop * 5))
        (drop, drop + path2.len() - path1.len())
    };

    for i in 1..path1.len() {
        // fill in first column of this row
        let this_cell_score = if path1[i] == path2[0] {
            *segment_lengths
                .get(&path1[i].abs())
                .ok_or(InversionError::SegmentNotFound(path1[i]))?
        } else {
            -1 * *segment_lengths
                .get(&path1[i].abs())
                .ok_or(InversionError::SegmentNotFound(path1[i]))?
        };
        let possible_scores = [0, -1, score_row_previous[0], -1];
        score_row_current[0] = amax(&possible_scores) + this_cell_score;
        let traceback_value = argmax(&possible_scores).try_into().expect("i32 overflow");
        if traceback_value != 0 {
            traceback_matrix.insert((i.try_into().expect("i32 overflow"), 0), traceback_value);
        }

        // fill in the rest of this row
        let len_i = *segment_lengths
            .get(&path1[i].abs())
            .ok_or(InversionError::SegmentNotFound(path1[i]))?;
        for j in 1..path2.len() {
            let len_j = *segment_lengths
                .get(&path2[j].abs())
                .ok_or(InversionError::SegmentNotFound(path2[j]))?;

            // heuristic: if we are too far from diagonal, leave traceback as implicit 0 and
            // calculate score as if we are starting alignment here regardless of what is in
            // cells nearby
            if ((i > j) && (i - j > max_row_drop)) || ((j > i) && (j - i > max_col_drop)) {
                score_row_current[j] = if path1[i] == path2[j] {
                    len_i
                } else {
                    -len_i - len_j
                };
            } else {
                let possible_scores = if path1[i] == path2[j] {
                    [
                        len_i,
                        score_row_previous[j - 1] + len_i, // come from diagonal
                        score_row_previous[j] + len_i,     // come from above
                        score_row_current[j - 1] + len_i,  // come from left
                    ]
                } else {
                    [
                        -len_i - len_j,
                        score_row_previous[j - 1] - len_i - len_j,
                        score_row_previous[j] - len_i,
                        score_row_current[j - 1] - len_j,
                    ]
                };

                score_row_current[j] = amax(&possible_scores);
                let traceback_value = argmax(&possible_scores).try_into().expect("i32 overflow");
                if traceback_value != 0 {
                    traceback_matrix.insert(
                        (
                            i.try_into().expect("i32 overflow"),
                            j.try_into().expect("i32 overflow"),
                        ),
                        traceback_value,
                    );
                }
            }
        }
        // update max/argmax of score matrix
        let (row_max, row_argmax) = max_and_argmax(&score_row_current);
        if row_max > max_score {
            max_score = row_max;
            argmax_score = (i.try_into().expect("i32 overflow"), row_argmax);
        }
        // now, switch rows. Not beautiful but faster than reallocating memory
        let score_row_tmp = score_row_previous;
        score_row_previous = score_row_current;
        score_row_current = score_row_tmp;
    }
    let traceback = traceback_lowmem(&path1, &path2, argmax_score, &traceback_matrix);
    debug!(
        "Finished lowmem alignment of length {}x{}",
        traceback.alignment_path1.len(),
        traceback.alignment_path2.len()
    );
    return Ok(traceback);
}

fn traceback_lowmem(
    path1: &[i32],
    path2: &[i32],
    argmax_score: (i32, i32),
    traceback_matrix: &HashMap<(i32, i32), i8>,
) -> align::Alignment {
    let (mut i, mut j) = argmax_score;
    let path1_end_index = i;
    let mut alignment_end_reached = false;
    let mut alignment_path1: Vec<i32> = Vec::new();
    let mut alignment_path2: Vec<i32> = Vec::new();
    while !alignment_end_reached {
        let segment_path1 = path1[<i32 as TryInto<usize>>::try_into(i).unwrap()];
        let segment_path2 = path2[<i32 as TryInto<usize>>::try_into(j).unwrap()];
        if alignment_path1.is_empty() || *alignment_path1.last().unwrap() != segment_path1 {
            alignment_path1.push(segment_path1);
        }
        if alignment_path2.is_empty() || *alignment_path2.last().unwrap() != segment_path2 {
            alignment_path2.push(segment_path2);
        }

        match traceback_matrix.get(&(i, j)) {
            Some(1) => {
                i -= 1;
                j -= 1;
            }
            Some(2) => i -= 1,
            Some(3) => j -= 1,
            None => alignment_end_reached = true,
            Some(x) => panic!("Bad value {} in traceback matrix!", x),
        }
    }
    let path1_start_index = i;
    alignment_path1.reverse();
    alignment_path2.reverse();
    return align::Alignment {
        alignment_path1,
        alignment_path2,
        path1_start_index,
        path1_end_index,
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_align_subproblem_lowmem() {
        let path1 = vec![2, 3, 4, -5, 6];
        let path2 = vec![6, 2, 7, -5];
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
        let alignment =
            align_paths_subproblem_lowmem(&path1, &path2, &segment_lengths, 100).unwrap();
        assert_eq!(alignment.alignment_path1, vec![2, 3, 4, -5]);
        assert_eq!(alignment.alignment_path2, vec![2, 7, -5]);
        assert_eq!(alignment.path1_start_index, 0);
        assert_eq!(alignment.path1_end_index, 3);
    }
}
