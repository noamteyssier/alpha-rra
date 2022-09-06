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
        Self::from_buffer(reader)
    }

    pub fn from_buffer<R: BufRead>(reader: R) -> Result<Self> {
        let (genes, pvalues) = reader
            .lines()
            .map(|line| line.expect("Malformed File, Unexpected EOL"))
            .enumerate()
            .filter_map(|(idx, line)| {
                match Self::parse_table(&line) {
                    Ok(r) => Some(r),
                    Err(e) => {
                        if idx == 0 {
                            None
                        } else {
                            panic!("Unable to parse input table: {}", e);
                        }
                    }
                }
            })
            .unzip();
        Ok(Self::from_raw_parts(genes, pvalues))
    }

    fn from_raw_parts(genes: Vec<String>, pvalues: Vec<f64>) -> Self {
        Self {
            genes,
            pvalues: Array1::from_vec(pvalues)
        }
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

#[cfg(test)]
mod testing {
    use std::io::Cursor;
    use super::Input;

    fn get_genes() -> Vec<String> {
        (0..10)
            .map(|x| format!("gene.{:?}", x))
            .collect()
    }
    
    fn get_pvalues() -> Vec<f64> {
        (1..=10)
            .map(|x| x as f64 / 1000.)
            .collect()
    }
    
    fn get_pvalues_string() -> Vec<String> {
        (1..=10)
            .map(|x| x as f64 / 1000.)
            .map(|x| format!("{}", x))
            .collect()
    }

    fn get_pvalues_string_malformed() -> Vec<String> {
        (1..=10)
            .map(|x| {
                if x != 5 { 
                    format!("{}", x as f64 / 1000.)
                } else {
                    String::from("NOT A FLOAT")
                }
            })
            .collect()
    }

    fn get_string_table(with_header: bool, malformed_pvalue: bool) -> String {
        let mut genes = get_genes();
        let mut pvalues = if malformed_pvalue {
            get_pvalues_string_malformed()
        } else {
            get_pvalues_string()
        };
        if with_header {
            genes.insert(0, "gene".to_string());
            pvalues.insert(0, "pvalues".to_string());
        }
        genes.iter().zip(pvalues.iter())
            .map(|(g, p)| format!("{}\t{}\n", g, p))
            .collect()
    }

    fn get_string_table_single_column() -> String {
        let genes = get_genes();
        genes.join("\n")
    }

    fn get_buffer_no_header() -> Cursor<String> {
        let string_table = get_string_table(false, false);
        Cursor::new(string_table)
    }

    fn get_buffer_with_header() -> Cursor<String> {
        let string_table = get_string_table(true, false);
        Cursor::new(string_table)
    }

    fn get_buffer_no_header_malformed() -> Cursor<String> {
        let string_table = get_string_table(false, true);
        Cursor::new(string_table)
    }

    fn get_buffer_single_column() -> Cursor<String> {
        Cursor::new(get_string_table_single_column())
    }

    #[test]
    fn test_from_raw_parts() {
        let genes = get_genes();
        let pvalues = get_pvalues();
        let input = Input::from_raw_parts(genes.clone(), pvalues.clone());
        assert_eq!(input.pvalues().to_vec(), pvalues);
        assert_eq!(input.genes().len(), genes.len());
    }

    #[test]
    fn test_from_buffer_no_header() {
        let buffer = get_buffer_no_header();
        let input = Input::from_buffer(buffer).unwrap();
        assert_eq!(input.pvalues().len(), 10);
        assert_eq!(input.genes().len(), 10);
    }

    #[test]
    fn test_from_buffer_with_header() {
        let buffer = get_buffer_with_header();
        let input = Input::from_buffer(buffer).unwrap();
        assert_eq!(input.pvalues().len(), 10);
        assert_eq!(input.genes().len(), 10);
    }

    #[test]
    #[should_panic]
    fn test_from_buffer_no_header_malformed_pvalue() {
        let buffer = get_buffer_no_header_malformed();
        Input::from_buffer(buffer).unwrap();
    }

    #[test]
    #[should_panic]
    fn test_from_buffer_no_header_single_column() {
        let buffer = get_buffer_single_column();
        Input::from_buffer(buffer).unwrap();
    }
}
