use structopt::StructOpt;
use std::path::PathBuf;

/// Efficient parrellel computation of the single cell Jaccard index
#[derive(StructOpt, Debug, Clone)]
#[structopt(name = "scJaccard")]
pub struct Opt {

    /// Read counts per region in MatrixMarket format. Any lines beginning with % will be treated
    /// as metadata. The first line of data is assumed to contain the totals for each column, and
    /// is used to preallocate memory. 
    #[structopt(short, long)]
    pub input: PathBuf,

    /// Regions of interest in the same order as specified in the MatrixMarket input in tab
    /// delimited bed format. No validation is done on chromosome or contig names, but they must be
    /// in the same format as used for the known bulk ATAC regions (chr1 versus 1, etc).
    #[structopt(short, long, parse(from_os_str))]
    pub bed: PathBuf,

    /// Cellular barcodes (or generalised identifiers) for each of the experiments listed in the
    /// MatrixMarket input. 
    #[structopt(long, parse(from_os_str))]
    pub barcodes: PathBuf,

    /// Sorted, tabix indexed bed files compressed with bgzip containing regions of interest for
    /// particular, potentially purified, bulk populations of cells. See the documentation for a
    /// discussion on how this data must be formatted, or use preformatted data distributed with
    /// avocato.
    #[structopt(short, long, parse(from_os_str))]
    pub atac: PathBuf,

    /*/// Output file, stdout if not present. 
    #[structopt(short, long, parse(from_os_str))]
    pub output: PathBuf,*/

    /// Number of threads to parrelise the computations over.  
    #[structopt(long, default_value = "1")]
    pub cores: usize,

    /// The number of valid contigs in given bed file. Useful if you are using an atypical genome build, or small
    /// test data sets.
    #[structopt(long, default_value = "22")]
    pub nchr: u64,

    /// Log level. Defaults to Info (useful information and statistics). 
    #[structopt(long, default_value = "info")]
    pub loglevel: String,

    /// Report full output (includes total number of intersections, isec, union, as well as
    /// jaccard). 
    #[structopt(long)]
    pub full: bool
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


/// Holds meta data for a particular cell.
pub struct Cell {
    /// Intersection
    pub isec: u32,
    pub union: u32,
    pub nisec: u32,
    pub nreg: u32
}

impl Cell { 
    pub fn jaccard(&self, known: f64) -> f64 { 
        let jac = self.isec as f64 / (self.union as f64 + known - self.isec as f64);
        return jac as f64
    }

}
