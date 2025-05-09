use regex::Regex;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

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
/// use inversion_finder::gfa;
///
/// assert_eq!(gfa::parse_gfa_path("1-,2+,3-,235+").unwrap(), vec![-1, 2, -3, 235]);
/// ```
pub fn parse_gfa_path(path_string: &str) -> Result<Vec<i32>, Box<dyn Error>> {
    let re = Regex::new(r"(\d+)([+-])").unwrap();
    let mut path_list: Vec<i32> = Vec::new();
    for segment in path_string.split(",") {
        let caps = re
            .captures(segment)
            .ok_or(format!("Invalid segment format '{}'", segment))?;
        path_list.push(
            match caps
                .get(2)
                .ok_or(format!("Invalid segment format '{}'", segment))?
                .as_str()
            {
                "+" => caps
                    .get(1)
                    .ok_or(format!("Invalid segment format '{}'", segment))?
                    .as_str()
                    .parse::<i32>()?,
                "-" => {
                    -1 * caps
                        .get(1)
                        .ok_or(format!("Invalid segment format '{}'", segment))?
                        .as_str()
                        .parse::<i32>()?
                }
                _ => return Err(format!("Invalid segment format '{}'", segment).into()),
            },
        )
    }
    return Ok(path_list);
}

/// Read a GFA into memory.
///
/// Only keep the information in the GFA that we will need later: the length of each segment, and
/// the paths.
///
/// # Arguments
/// * `path` - path to GFA to parse
///
/// # Returns
/// * `segment_lengths` - map of segment ID to segment length in bp
/// * `paths` - map of path ID to vector of path with segment orientation indicated by sign
/// * `path_names` - the keys of `paths`, but in the order they were read, because I find the
///   nondeterministic order that stuff comes out of the HashMap to be disturbing
pub fn read_gfa(
    path: PathBuf,
) -> Result<(HashMap<i32, i32>, HashMap<String, Vec<i32>>, Vec<String>), Box<dyn Error>> {
    let file = match File::open(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.display(), why),
        Ok(file) => file,
    };
    let reader = BufReader::new(file);

    let mut segment_lengths: HashMap<i32, i32> = HashMap::new();
    let mut paths = HashMap::new();
    let mut path_names = Vec::new();

    for line_result in reader.lines() {
        let line = line_result?;
        let fields: Vec<&str> = line.split("\t").collect();

        if fields[0] == "S" {
            segment_lengths.insert(fields[1].to_string().parse()?, fields[2].len().try_into()?);
        } else if fields[0] == "P" {
            paths.insert(fields[1].to_string(), parse_gfa_path(fields[2])?);
            path_names.push(fields[1].to_string());
        }
    }

    return Ok((segment_lengths, paths, path_names));
}

/// Lookup start and end positions of segments in a path.
///
/// # Arguments
///
/// * `path`: list of segments in path
/// * `segment_lengths`: map of segment ID to segment length in bp
/// * `segments`: indices in `path` of segments to look up
///
/// # Returns
///
/// * map of segment index to tuple of start and end positions of segment in path, in bp
///
/// # Examples
///
/// ```
/// use inversion_finder::gfa;
/// use std::collections::HashMap;
///
/// let mut segment_lengths: HashMap<i32, i32> = HashMap::new();
/// for i in 0..8 {
///     segment_lengths.insert(i, 100);
/// }
/// let base_positions =
///     gfa::lookup_base_positions(&[1,-2,3,4,5,-6,-7], &segment_lengths, &[1,3,6]).unwrap();
/// assert_eq!(base_positions.get(&1).unwrap(), &(101, 200));
/// assert_eq!(base_positions.get(&6).unwrap(), &(601, 700));
/// ```
pub fn lookup_base_positions(
    path: &[i32],
    segment_lengths: &HashMap<i32, i32>,
    segment_indices: &[i32],
) -> Result<HashMap<i32, (i32, i32)>, Box<dyn Error>> {
    let segment_indices_set = HashSet::<i32>::from_iter(segment_indices.iter().cloned());
    let mut segment_positions: HashMap<i32, (i32, i32)> = HashMap::new();

    let mut current_position = 0;
    for (i, segment) in path.iter().map(|s| s.abs()).enumerate() {
        let this_segment_length = segment_lengths
            .get(&segment)
            .ok_or(format!("Cannot find segment {}", segment))?;
        if segment_indices_set.contains(&i.try_into()?) {
            segment_positions.insert(
                i.try_into()?,
                (current_position + 1, current_position + this_segment_length),
            );
        }
        current_position += this_segment_length;
    }
    return Ok(segment_positions);
}

mod tests {
    #[test]
    #[should_panic]
    fn test_parse_bad_gfa_path() {
        super::parse_gfa_path("1-,2+,3-,235").unwrap();
    }
}
