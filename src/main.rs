//! `scJaccard` is a tool for quantifying the overlap between single cell peak calls and known
//! purified cell types. 
//!
//! To use this tool, you need:
//!     - Peak calls in triplet MatrixMarket format (see
//!     [here](https://math.nist.gov/MatrixMarket/formats.html) for details)
//!     - Known peaks in `.bed` format, sorted and index by `tabix`.
//!
//! To run the tool, simply give the paths as arguments.
//!
//! ```sh
//! scJaccard --input   path/to/matrixMarket \
//!           --bed     path/to/file.bed \
//!           --atac    path/to/atac.bed
//! ```
//!
//! ## Input Formats
//!
//! 1. MatrixMarket. 
//! 2. Indexed `bed`.
//!
//! ## Associated Tools
//!
//! This tool has been created as a part of the `avocato` single cell Assay for
//! Transposase-Accessible Chromatin using sequencing (scATAC-seq) pipeline.  
//!

// This is a bit overkill because it turns off
// linting for the whole program, but I hate this
// warning.
#![allow(non_snake_case)]

//use std::io::Write;
use structopt::StructOpt;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::collections::HashMap;
use log::info;
use env_logger::Env;
use rust_htslib::tbx::{self, Read};
use rayon::prelude::*;
use std::fmt::Write as OtherWrite;

mod utils;
mod structs;

fn main() -> std::io::Result<()> {
    // Parse options using structopt
    // Requires update of crates.
    let opt = structs::Opt::from_args();

    // I want a logger, but I want to 
    // control it.
    let env = Env::default()
        .filter_or("MY_LOG_LEVEL", opt.loglevel)
        .write_style_or("MY_LOG_STYLE", "always");
    env_logger::init_from_env(env);

    let mut stats = HashMap::new();

    // Parse the barcodes tsv file.
    info!("Reading barcodes.");
    let mut barcodes = HashMap::new();
    {
       let bcodes = File::open(&opt.barcodes)?;
       let bcodes_reader = BufReader::new(bcodes);

       let mut r = 1;
       for line in bcodes_reader.lines() {
           barcodes.insert(r.to_string(),
            line?);
           r += 1
       }
    }

    // Parse the regions in the bed file
    // to create a hashmap of regions 
    // (chr, start, stop) indexed by
    // their order, as is referred to in 
    // the Matrix file.
    info!("Reading bed file coordinates.");
    let mut regions = HashMap::new();
    {
        let bed_file = File::open(&opt.bed)?;
        let bed_file_reader = BufReader::new(bed_file);

        let mut r = 1;
        for line in bed_file_reader.lines() {
            // Parse the line.
            let _line = line?;
            let vec = _line.split("\t").collect::<Vec<&str>>();
            // Insert a new record with the numeric index as 
            // the hashed portion.
            regions.insert(r.to_string(), 
                       structs::Region {
                           chr: vec[0].to_string(),
                           start: vec[1].parse::<u32>().unwrap(),
                           stop: vec[2].parse::<u32>().unwrap() });
            r += 1
        }

    }

    // Figure out how many cells there are
    // in order to preallocate that memory
    // for the HashMap.
    let n_cells: usize = utils::count_cells(&opt.input)?;

    // Parse the MatrixMarket format file
    // to obtain a hashmap of region (indices)
    // indexed by cell barcode (indices)A
    //
    // Open an empty scope to create,
    // effectively, a context handler
    // for the file connection.
    
    info!("Reading MatrixMarket data.");
    let mut first: bool = true; 
    let mut mtx = HashMap::with_capacity(n_cells);
    {
        let file = File::open(&opt.input)?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let vec = line?;
            if !vec.starts_with('%') {
                // Skip the first numeric line as we 
                // use it for pre-allocation.
                if first { 
                    first = false;
                    continue; 
                } 
                let vec = vec.split(" ");
                let vec = vec.collect::<Vec<&str>>();
                if !&mtx.contains_key(&vec[1].to_string()) {
                    let idx = vec[1].to_string();
                    mtx.insert(vec[1].to_string(), 
                               Vec::<structs::Record>::new());

                    mtx.get_mut(&idx).unwrap().push(structs::Record {
                                    i: vec[0].parse::<u32>().unwrap(), 
                                    //j: vec[1].parse::<u32>().unwrap(), 
                                    //val: vec[2].parse::<u32>().unwrap()
                                    });
                    stats.insert(vec[1].to_string(),
                                 structs::Cell { 
                                     isec: 0, 
                                     union: regions[&vec[0].to_string()].total() });
                                        
                } else { 
                    let idx = vec[1].to_string();
                    mtx
                        .get_mut(&idx)
                        .unwrap()
                        .push(structs::Record {
                            i: vec[0].parse::<u32>().unwrap()
                            //j: vec[1].parse::<u32>().unwrap(), 
                            //val: vec[2].parse::<u32>().unwrap()
                            });
                    stats
                        .get_mut(&vec[1].to_string())
                        .unwrap()
                        .union += regions[&vec[0].to_string()].total();
                }

                }

            } 
    } 
    //println!("Number of records for 1: {}", stats["1"].union);
    info!("Reading known regions.");
    let path_bed = opt.atac;
    let mut tbx_reader = tbx::Reader::from_path(&path_bed).unwrap();

    // Holder for the totals.
    let mut known = structs::Cell{ 
        isec: 0, 
        union: 0 
    };

    let nchr: u64 = opt.nchr;

    // Go through the records once and get the sum 
    // There is no .fetch_string() impl for tbx::Reader
    // so I'm just going to go through the chromosomes manually and
    // do it.
    for i in 0..nchr { 
        tbx_reader
            .fetch(i, 0, 1000000000).unwrap(); // All of it?
        for record in tbx_reader.records() { 
            let r = utils::parse_bed_line(record.unwrap());
            known
                .union += r.total();
        }
    }

    //println!("total {}", known.union);

    info!("Using {} threads.", rayon::current_num_threads());
    let isec: HashMap<String, u32> = mtx.par_iter()
        .map(|(key,value)| utils::total_isec(key, value, &path_bed, &regions))
        .collect::<HashMap<String, u32>>();

    for (key, value) in isec {
        stats
            .get_mut(&key)
            .unwrap()
            .isec += value
    }

    info!("Computed jaccard for {} cells.", stats.len());
    //info!("Writing output to {:?}.", &opt.output);

    //let mut output = File::create(opt.output)?; 
    for (key, value) in stats {
        let mut s = String::new();
        write!(s, "{} {:.5}", barcodes[&key.to_string()], value.jaccard(*&known.union as f64)).unwrap();
        println!("{}", s);
        //output.write(s.as_bytes())?;
    }

    info!("scJaccard finished successfully.");
    Ok(())
}

