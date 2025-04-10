pub mod align;
pub mod gfa;
pub mod lowmem;

/// Find the maximum value in a list.
///
/// # Examples
///
/// ```
/// assert_eq!(inversion_finder::max(&[2,5,3]), 5);
/// ```
pub fn max(a: &[i32]) -> i32 {
    *a.iter().max().unwrap()
}

/// Find the index of the maximum value in a list.
///
/// # Examples
///
/// ```
/// assert_eq!(inversion_finder::argmax(&[2,5,3]), 1);
/// ```
pub fn argmax(a: &[i32]) -> i32 {
    (0..a.len())
        .max_by_key(|x| a[*x])
        .unwrap()
        .try_into()
        .unwrap()
}
