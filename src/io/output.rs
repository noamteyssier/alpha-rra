use std::{io::{Write, stdout}, fs::File};
use anyhow::Result;
use alpha_rra::ResultsRRA;

/// Match the output stream on a path name
pub fn match_writer(path: &str) -> Box<dyn Write> {
    if let Ok(file) = File::create(path) {
        Box::new(file)
    } else {
        Box::new(stdout())
    }
}

/// Writes the results of the analysis to a tab-delim `Write` object
pub fn write_lines(writer: &mut Box<dyn Write>, results: &ResultsRRA) -> Result<()> {
    writeln!(writer, "{}\t{}\t{}\t{}", "name", "score", "pvalue", "adj_pvalue")?;
    for (name, score, pvalue, adj_pvalue) in results.zip() {
        writeln!(writer, "{}\t{}\t{}\t{}", name, score, pvalue, adj_pvalue)?
    }
    Ok(())
}
