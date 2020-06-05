use std::path::PathBuf;



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
pub fn count_cells (path: &PathBuf) -> usize {
    let file = File::open(opt.input)?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let vec = line?;
        if !vec.starts_with('%') {
            vec[1].parse<usize>()

        }

}
