use std::path::PathBuf;
use structopt::StructOpt;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::collections::HashMap;

/// A basic example
#[derive(StructOpt, Debug)]
#[structopt(name = "basic")]
struct Opt {

    /// Input file
    #[structopt(short, long)]
    input: PathBuf,

    /// Output file
    #[structopt(short, long, parse(from_os_str))]
    output: PathBuf,
}

struct Region {
    i: u32,
    j: u32,
    val: u32
}


fn main() -> std::io::Result<()> {
    let opt = Opt::from_args();

    let mut regions = HashMap::new();

    let file = File::open(opt.input)?;
    let reader = BufReader::new(file);

    let mut i = 0;

    for line in reader.lines() {
        //regions.insert(line[0], line[1]);
        //

        let vec = line?;
        if !vec.starts_with('%') {
            let vec = vec.split(" ");
            let vec = vec.collect::<Vec<&str>>();
            //println!("{}", vec[0]);

            if !&regions.contains_key(&vec[2].to_string()) {
                regions.insert(vec[2].to_string(), Region {i: vec[0].parse::<u32>().unwrap(), 
                    j: vec[1].parse::<u32>().unwrap(), 
                    val: vec[2].parse::<u32>().unwrap()});
            }

            for (key, value) in &regions {
                println!("{} / {}", key, value.j);
            }
        }

        i = i + 1;

    }

    Ok(())
}
