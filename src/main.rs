use std::path::PathBuf;
use structopt::StructOpt;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::collections::HashMap;
use log::info;
use env_logger::Env;
use rust_htslib::tbx::{self, Read};
use rayon::prelude::*;
use std::sync::Arc;

#[derive(StructOpt, Debug)]
#[structopt(name = "basic")]
struct Opt {

    /// Input file
    #[structopt(short, long)]
    input: PathBuf,

    /// Output file
    #[structopt(short, long, parse(from_os_str))]
    output: PathBuf,

    /// Bed file
    #[structopt(short, long, parse(from_os_str))]
    bed: PathBuf,

    /// Known ATAC peaks
    #[structopt(short, long, parse(from_os_str))]
    atac: PathBuf
}

struct Record {
    i: u32
    //j: u32,
    //val: u32
}

#[derive(Clone)]
struct Region {
    chr: String,
    start: u32,
    stop: u32
}

impl Region { 
    fn total(&self) -> u32 { 
        return self.stop - self.start
    }
}

struct Cell {
    isec: u32,
    union: u32
}

fn get_intersection(interval1: &Region, 
                    interval2: &Region) -> u32 {

    if interval2.start > interval1.stop || interval1.start > interval2.stop {
        return 0; // does not contribute
    }

    let start = std::cmp::max(interval1.start, interval2.start);
    let stop = std::cmp::min(interval1.stop, interval2.stop);

    return stop - start;
}

fn parse_bed_line(line: Vec<u8>) -> Region {
    let r = String::from_utf8(line).unwrap();
    let r = r.split("\t").collect::<Vec<&str>>();

    return Region {
        chr: r[0].to_string(),
        start: r[1].parse::<u32>().unwrap(),
        stop: r[2].parse::<u32>().unwrap()
    }
}

fn total_isec(cell_barcode: &String, 
              records: &Vec<Record>,
              path_bed: &PathBuf,
              regions: &HashMap<String, Region>) -> (String, u32) {

    let mut tbx_reader = tbx::Reader::from_path(path_bed).unwrap();
    let mut isec = 0;
    for r in records { 
        let reg = &regions[&r.i.to_string()];
        let tid = tbx_reader.tid(&reg.chr).unwrap();
        tbx_reader.fetch(tid, reg.start as u64, reg.stop as u64).unwrap();


        for record in tbx_reader.records() {
            let parsed_region = parse_bed_line(record.unwrap());
            isec += get_intersection(reg, &parsed_region);
        }

    }

    let bar = cell_barcode.clone();

    return (bar, isec as u32);
    
}

fn main() -> std::io::Result<()> {
    // Parse options using structopt
    // Requires update of crates.
    let opt = Opt::from_args();

    // I want a logger, but I want to 
    // control it.
    let env = Env::default()
        .filter_or("MY_LOG_LEVEL", "info")
        .write_style_or("MY_LOG_STYLE", "always");
    env_logger::init_from_env(env);

    let mut stats = HashMap::new();

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
                       Region {
                           chr: vec[0].to_string(),
                           start: vec[1].parse::<u32>().unwrap(),
                           stop: vec[2].parse::<u32>().unwrap() });
            r += 1
        }

    }

    // Parse the MatrixMarket format file
    // to obtain a hashmap of region (indices)
    // indexed by cell barcode (indices)A
    //
    // Open an empty scope to create,
    // effectively, a context handler
    // for the file connection.
    
    info!("Reading MatrixMarket data.");
    let mut mtx = HashMap::new();
    {
        let file = File::open(opt.input)?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let vec = line?;
            //println!("{:?}", vec);
            if !vec.starts_with('%') {
                let vec = vec.split(" ");
                let vec = vec.collect::<Vec<&str>>();
                if !&mtx.contains_key(&vec[1].to_string()) {
                    let idx = vec[1].to_string();
                    mtx.insert(vec[1].to_string(), 
                               Vec::<Record>::new());

                    mtx.get_mut(&idx).unwrap().push(Record {
                                    i: vec[0].parse::<u32>().unwrap(), 
                                    //j: vec[1].parse::<u32>().unwrap(), 
                                    //val: vec[2].parse::<u32>().unwrap()
                                    });
                    stats.insert(vec[1].to_string(),
                                 Cell { 
                                     isec: 0, 
                                     union: regions[&vec[0].to_string()].total() });
                                        
                } else { 
                    let idx = vec[1].to_string();
                    mtx
                        .get_mut(&idx)
                        .unwrap()
                        .push(Record {
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
    println!("Number of records for 1: {}", stats["1"].union);
    info!("Reading known regions.");
    let path_bed = opt.atac;
    let mut tbx_reader = tbx::Reader::from_path(&path_bed).unwrap();

    // Holder for the totals.
    let mut known = Cell{ 
        isec: 0, 
        union: 0 
    };

    // Go through the records once and get the sum 
    // There is no .fetch_string() impl for tbx::Reader
    // so I'm just going to go through the chromosomes manually and
    // do it.
    for i in 1..22 { 
        tbx_reader
            .fetch(i, 0, 1000000000).unwrap(); // All of it?
        for record in tbx_reader.records() { 
            let r = parse_bed_line(record.unwrap());
            known
                .union += r.total();
        }
    }

    println!("total {}", known.union);

    println!("Threads: {}", rayon::current_num_threads());
    let isec: HashMap<String, u32> = mtx.par_iter()
        .map(|(key,value)| total_isec(key, value, &path_bed, &regions))
        .collect::<HashMap<String, u32>>();


    /*for (key, value) in mtx { 
        for r in value{
            let reg = &regions[&r.i.to_string()];
            let tid = tbx_reader.tid(&reg.chr).unwrap();
            tbx_reader.fetch(tid, reg.start as u64, reg.stop as u64).unwrap();

            for record in tbx_reader.records() {
                let parsed_region = parse_bed_line(record.unwrap());
                //let reg = reg;
                stats
                    .get_mut(&key)
                    .unwrap()
                    .isec += get_intersection(reg, &parsed_region);
            }

    }}
    */

    for (key, value) in stats {
        println!("Isec: {}, Union: {}", value.isec, value.union);
    }

    Ok(())
}

