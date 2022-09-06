use alpha_rra::alpha_rra;
use clap::Parser;
use anyhow::Result;
use adjustp::Procedure;

mod io;
use io::{Input, match_writer};

use crate::io::write_lines;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {

    /// Filepath of the input count matrix
    #[clap(short, long, value_parser)]
    input: String,

    /// Permutations
    #[clap(short, long, value_parser, default_value="100")]
    permutations: usize,

    /// Alpha Threshold
    #[clap(short, long, value_parser, default_value="0.1")]
    alpha: f64,

    /// Multiple Hypothesis Correction (bonferroni, bh, by)
    #[clap(short='f', long, value_parser, default_value="bh")]
    correction: String
}

fn main() -> Result<()> {
    let args = Args::parse();
    let input = Input::from_path(&args.input)?;
    let correction = match args.correction.as_str() {
        "bonferroni" => Procedure::Bonferroni,
        "bh" | "fdr" => Procedure::BenjaminiHochberg,
        "by" => Procedure::BenjaminiYekutieli,
        _ => panic!("Unexpected correction method provided: {}", args.correction)
    };

    let results = alpha_rra(
        input.pvalues(), 
        input.genes(), 
        args.alpha, 
        args.permutations,
        correction
    );

    let mut output = Box::new(stdout());
    write_lines(&mut output, &results)?;

    Ok(())
}
