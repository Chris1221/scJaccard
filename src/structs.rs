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

   /* /// Output file, stdout if not present. 
    #[structopt(short, long, parse(from_os_str))]
    pub output: PathBuf, */

    /// Number of cores to use. Defaults to 1. 
    #[structopt(long, default_value = "1")]
    pub cores: usize,

    /// Number of contigs in the ATAC file to iterate over. Defaults to 22. 
    #[structopt(long, default_value = "22")]
    pub nchr: u64,

    /// Log level. Defaults to Info (useful information and statistics). 
    #[structopt(long, default_value = "info")]
    pub loglevel: String
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
