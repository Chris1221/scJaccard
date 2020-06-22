use std::io::{BufRead, BufReader};
use std::fs::File;
use std::path::PathBuf;
use std::collections::HashMap;
use rust_htslib::tbx::{self, Read};
use crate::structs;
use log::debug;

/// Determines the number of cells to pre-allocate space 
/// for a HashMap.
///
/// Reads from the given buffer the first non-commented
/// line and returns the total number of j elements
/// as a `usize`.
///
/// # Examples
///
/// Given a `path: &PathBuf`.
///
/// ```
/// let n: usize = count_cells(path);
/// ```
pub fn count_cells (path: &PathBuf) -> Result<usize, std::io::Error> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut n: usize = 0;

    for line in reader.lines() {
        let vec = line?;
        if !vec.starts_with('%') {
            let vec = vec.split(" ");
            let vec = vec.collect::<Vec<&str>>();
            n = vec[1].parse::<usize>().unwrap();
        }
    }
    Ok(n)
}

pub fn get_intersection(interval1: &structs::Region, 
                    interval2: &structs::Region) -> u32 {

    if interval2.start > interval1.stop || interval1.start > interval2.stop {
        return 0; // does not contribute
    }

    let start = std::cmp::max(interval1.start, interval2.start);
    let stop = std::cmp::min(interval1.stop, interval2.stop);

    return stop - start;
}

pub fn parse_bed_line(line: Vec<u8>) -> structs::Region {
    let r = String::from_utf8(line).unwrap();
    let r = r.split("\t").collect::<Vec<&str>>();

    return structs::Region {
        chr: r[0].to_string(),
        start: r[1].parse::<u32>().unwrap(),
        stop: r[2].parse::<u32>().unwrap()
    }
}

pub fn total_isec(cell_barcode: &String, 
              records: &Vec<structs::Record>,
              path_bed: &PathBuf,
              regions: &HashMap<String, structs::Region>) -> (String, (u32, u32, u32)) {

    let mut tbx_reader = tbx::Reader::from_path(path_bed).unwrap();
    let mut isec = 0;
    let mut nint = 0;
    let mut nreg = 0;
    for r in records { 
        let reg = &regions[&r.i.to_string()];
        let tid = tbx_reader.tid(&reg.chr).unwrap();
        let s: u64;
        if reg.start > 10000 {
            s = reg.start as u64 - 10000;
        } else {
            s = reg.start as u64;
        }

        tbx_reader.fetch(tid, s, reg.stop as u64 + 10000).unwrap();

        debug!("* Region start, stop ({}, {})", reg.start, reg.stop); 

        for record in tbx_reader.records() {
            let parsed_region = parse_bed_line(record.unwrap());
            debug!("** Returned regions {}, {}", parsed_region.start, parsed_region.stop);
            let isecs = get_intersection(reg, &parsed_region);
            if isecs > 0 {
                nint += 1
            }
            isec += isecs; 
            //debug!("*** ISEC now {}", isec);
            //nint += 1;
        }
        nreg += 1

    }

    let bar = cell_barcode.clone();

    return (bar, (isec as u32, nint as u32, nreg as u32) );
    
}
