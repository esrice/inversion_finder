use std::{error, fmt};

pub mod align;
pub mod alignment_interface;
pub mod gfa;
pub mod lowmem;

/// Find the maximum value in a list.
///
/// # Examples
///
/// ```
/// assert_eq!(inversion_finder::amax(&[2,5,3]), 5);
/// ```
pub fn amax(a: &[i32]) -> i32 {
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

#[derive(Debug)]
pub enum InversionError {
    GfaParse(String),
}

impl fmt::Display for InversionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            InversionError::GfaParse(e) => write!(f, "Error parsing GFA: {}", e),
        }
    }
}

impl error::Error for InversionError {}
