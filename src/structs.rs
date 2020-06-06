use structopt::StructOpt;
use std::path::PathBuf;

/// Input options from the Opt crate.
#[derive(StructOpt, Debug, Clone)]
#[structopt(name = "scJaccard_args")]
pub struct Opt {

    /// Input file
    #[structopt(short, long)]
    pub input: PathBuf,

    /// Bed file
    #[structopt(short, long, parse(from_os_str))]
    pub bed: PathBuf,

    /// Barcodes 
    #[structopt(long, parse(from_os_str))]
    pub barcodes: PathBuf,

    /// Known ATAC peaks
    #[structopt(short, long, parse(from_os_str))]
    pub atac: PathBuf,

    /// Output path 
    #[structopt(short, long, parse(from_os_str))]
    pub output: PathBuf
}

pub struct Record {
    pub i: u32
}

#[derive(Clone)]
pub struct Region {
    pub chr: String,
    pub start: u32,
    pub stop: u32
}

impl Region { 
    pub fn total(&self) -> u32 { 
        return self.stop - self.start
    }
}

pub struct Cell {
    pub isec: u32,
    pub union: u32
}

impl Cell { 
    pub fn jaccard(&self, known: f64) -> f64 { 
        let jac = self.isec as f64 / (self.union as f64 + known - self.isec as f64);
        return jac as f64
    }

}
