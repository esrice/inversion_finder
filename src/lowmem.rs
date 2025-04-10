use super::{amax, argmax};
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
    return (max, argmax.try_into().unwrap());
}

fn initialize_matrices_lowmem(
    path1: &[i32],
    path2: &[i32],
    segment_lengths: &HashMap<i32, i32>,
) -> (Vec<i32>, Vec<i32>, HashMap<(i32, i32), i8>, i32, (i32, i32)) {
    let mut score_row_previous = Vec::<i32>::with_capacity(path2.len());
    // this is a sparse representation of the traceback matrix: if a key isn't present, that
    // indicates that the value is actually 0
    let mut traceback_matrix = HashMap::<(i32, i32), i8>::new();

    // fill in the corner
    score_row_previous.push(if path1[0] == path2[0] {
        *segment_lengths.get(&path1[0].abs()).unwrap()
    } else {
        -1
    });

    // fill in the rest of the first row
    for j in 1..path2.len() {
        let this_cell_score = if path2[j] == path1[0] {
            *segment_lengths.get(&path2[j].abs()).unwrap()
        } else {
            -1
        };

        let possible_scores = [0, -1, -1, score_row_previous[j - 1]];
        score_row_previous.push(amax(&possible_scores) + this_cell_score);
        let traceback_value = argmax(&possible_scores).try_into().unwrap();
        if traceback_value != 0 {
            traceback_matrix.insert((0, j.try_into().unwrap()), traceback_value);
        }
    }

    let (max_score, argmax_score_j) = max_and_argmax(&score_row_previous);
    let argmax_score: (i32, i32) = (0, argmax_score_j);

    let mut score_row_current = Vec::<i32>::with_capacity(path2.len());
    for _ in 0..path2.len() {
        score_row_current.push(0);
    }

    return (
        score_row_previous,
        score_row_current,
        traceback_matrix,
        max_score,
        argmax_score,
    );
}

pub fn align_paths_subproblem_lowmem(
    path1: &[i32],
    path2: &[i32],
    segment_lengths: &HashMap<i32, i32>,
    drop: usize,
) -> (Vec<i32>, Vec<i32>) {
    let (
        mut score_row_previous,
        mut score_row_current,
        mut traceback_matrix,
        mut max_score,
        mut argmax_score,
    ) = initialize_matrices_lowmem(&path1, &path2, &segment_lengths);

    let (max_row_drop, max_col_drop) = if path1.len() > path2.len() {
        (drop + path1.len() - path2.len(), drop)
    } else {
        (drop, drop + path2.len() - path1.len())
    };

    for i in 1..path1.len() {
        // fill in first column of this row
        let this_cell_score = if path1[i] == path2[0] {
            *segment_lengths.get(&path1[i].abs()).unwrap()
        } else {
            -1
        };
        let possible_scores = [0, -1, score_row_previous[0], -1];
        score_row_current[0] = amax(&possible_scores) + this_cell_score;
        let traceback_value = argmax(&possible_scores).try_into().unwrap();
        if traceback_value != 0 {
            traceback_matrix.insert((i.try_into().unwrap(), 0), traceback_value);
        }

        // fill in the rest of this row
        for j in 1..path2.len() {
            let this_cell_score = if path1[i] == path2[j] {
                *segment_lengths.get(&path1[i].abs()).unwrap()
            } else {
                -1
            };

            // heuristic: if we are too far from diagonal, leave traceback as implicit 0 and
            // calculate score as if we are starting alignment here regardless of what is in
            // cells nearby
            if ((i > j) && (i - j > max_row_drop)) || ((j > i) && (j - i > max_col_drop)) {
                score_row_current[j] = this_cell_score;
            } else {
                let possible_scores = [
                    0,
                    score_row_previous[j - 1], // come from diagonal
                    score_row_previous[j],     // come from above
                    score_row_current[j - 1],  // come from left
                ];

                score_row_current[j] = amax(&possible_scores) + this_cell_score;
                let traceback_value = argmax(&possible_scores).try_into().unwrap();
                if traceback_value != 0 {
                    traceback_matrix.insert(
                        (i.try_into().unwrap(), j.try_into().unwrap()),
                        traceback_value,
                    );
                }
            }
        }
        // update max/argmax of score matrix
        let (row_max, row_argmax) = max_and_argmax(&score_row_current);
        if row_max > max_score {
            max_score = row_max;
            argmax_score = (i.try_into().unwrap(), row_argmax);
        }
        // now, switch rows. Not beautiful but faster than reallocating memory
        let score_row_tmp = score_row_previous;
        score_row_previous = score_row_current;
        score_row_current = score_row_tmp;
    }

    return traceback_lowmem(&path1, &path2, argmax_score, &traceback_matrix);
}

fn traceback_lowmem(
    path1: &[i32],
    path2: &[i32],
    argmax_score: (i32, i32),
    traceback_matrix: &HashMap<(i32, i32), i8>,
) -> (Vec<i32>, Vec<i32>) {
    let (mut i, mut j) = argmax_score;
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
    alignment_path1.reverse();
    alignment_path2.reverse();
    return (alignment_path1, alignment_path2);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_align_subproblem_lowmem() {
        let path1 = vec![2, 3, 4, -5, 6];
        let path2 = vec![6, 2, 7, -5];
        let mut segment_lengths: HashMap<i32, i32> = HashMap::new();
        for i in 0..7 {
            segment_lengths.insert(i, 100);
        }
        let (path1_alignment, path2_alignment) =
            align_paths_subproblem_lowmem(&path1, &path2, &segment_lengths, 100);
        assert_eq!(path1_alignment, vec![2, 3, 4, -5]);
        assert_eq!(path2_alignment, vec![2, 7, -5]);
    }
}
