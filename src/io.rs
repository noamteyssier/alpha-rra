use std::{fs::File, io::{BufReader, BufRead}, str::SplitWhitespace};
use anyhow::{Result, bail};
use ndarray::Array1;

pub struct Input {
    genes: Vec<String>,
    pvalues: Array1<f64>
}
impl Input {
    
    pub fn from_path(path: &str) -> Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let (genes, pvalues): (Vec<String>, Vec<f64>) = reader
            .lines()
            .map(|line| line.expect("Malformed File, Unexpected EOL"))
            .map(|line| Self::parse_table(&line).expect("Error processing row"))
            .unzip();
        Ok(Self {
            genes,
            pvalues: Array1::from_vec(pvalues)
        })
    }

    fn parse_table(row: &str) -> Result<(String, f64)> {
        let mut tokens = row.split_whitespace();
        Ok((
            Self::parse_gene(&mut tokens)?,
            Self::parse_pvalue(&mut tokens)?
        ))
    }

    fn parse_gene(tokens: &mut SplitWhitespace) -> Result<String> {
        if let Some(g) = tokens.next() {
            Ok(g.to_string())
        } else {
            bail!("Row is unexpectedly empty!")
        }
    }

    fn parse_pvalue(tokens: &mut SplitWhitespace) -> Result<f64> {
        if let Some(pvalue_string) = tokens.next() {
            if let Ok(pvalue) = pvalue_string.parse::<f64>() {
                Ok(pvalue)
            } else {
                bail!("Second column must be parseable into float!")
            }
        } else {
            bail!("Row must be at least 2 columns!")
        }
    }

    pub fn pvalues(&self) -> &Array1<f64> {
        &self.pvalues
    }

    pub fn genes(&self) -> &Vec<String> {
        &self.genes
    }

}
