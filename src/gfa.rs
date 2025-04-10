use regex::Regex;
use std::collections::HashMap;
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
/// assert_eq!(gfa::parse_gfa_path("1-,2+,3-,235+"), vec![-1, 2, -3, 235]);
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
///
/// # Arguments
/// * `path` - path to GFA to parse
///
/// # Returns
/// * `segment_lengths` - map of segment ID to segment length in bp
/// * `paths` - map of path ID to vector of path with segment orientation indicated by sign
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

/// Lookup start and end positions of segments in a path.
///
/// # Arguments
///
/// * `path`: list of segments in path
/// * `segment_lengths`: map of segment ID to segment length in bp
/// * `segments`: IDs of segments to look up
///
/// # Returns
///
/// * map of segment ID to tuple of start and end positions of segment in path, in bp
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
///     gfa::lookup_base_positions(&[1,-2,3,4,5,-6,-7], &segment_lengths, &[1,3,7]);
/// assert_eq!(base_positions.get(&1).unwrap(), &(1, 100));
/// assert_eq!(base_positions.get(&7).unwrap(), &(601, 700));
/// ```
pub fn lookup_base_positions(
    path: &[i32],
    segment_lengths: &HashMap<i32, i32>,
    segments: &[i32],
) -> HashMap<i32, (i32, i32)> {
    // first, populate the map of segment positions with the segments we actually want to lookup as
    // keys and dummy values
    let mut segment_positions: HashMap<i32, (i32, i32)> = HashMap::new();
    for segment in segments {
        segment_positions.insert(segment.abs(), (-1, -1));
    }

    let mut current_position = 0;
    for segment in path.iter().map(|s| s.abs()) {
        let this_segment_length = match segment_lengths.get(&segment) {
            Some(x) => x,
            None => panic!("{segment} not in segment_lengths!"),
        };
        if segment_positions.contains_key(&segment) {
            segment_positions.insert(
                segment,
                (current_position + 1, current_position + this_segment_length),
            );
        }
        current_position += this_segment_length;
    }
    return segment_positions;
}

mod tests {
    #[test]
    #[should_panic]
    fn test_parse_bad_gfa_path() {
        super::parse_gfa_path("1-,2+,3-,235");
    }
}
