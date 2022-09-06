use std::io::Write;
use anyhow::Result;
use alpha_rra::ResultsRRA;

/// Writes the results of the analysis to a tab-delim `Write` object
pub fn write_lines(writer: &mut Box<impl Write>, results: &ResultsRRA) -> Result<()> {
    writeln!(writer, "{}\t{}\t{}\t{}", "name", "score", "pvalue", "adj_pvalue")?;
    for (name, score, pvalue, adj_pvalue) in results.zip() {
        writeln!(writer, "{}\t{}\t{}\t{}", name, score, pvalue, adj_pvalue)?
    }
    Ok(())
}
