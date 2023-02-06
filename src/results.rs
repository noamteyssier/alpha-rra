use adjustp::{adjust, Procedure};

/// Type alias for the results of an Alpha-RRA run for a single entry
type EntryRRA<'a> = (&'a str, f64, f64, f64);

/// Handles the results of an Alpha-RRA run
pub struct ResultsRRA {
    names: Vec<String>,
    scores: Vec<f64>,
    pvalues: Vec<f64>,
    adj_pvalues: Vec<f64>,
}

impl ResultsRRA {
    /// Creates a new instance
    pub fn new(
        names: Vec<String>,
        scores: Vec<f64>,
        pvalues: Vec<f64>,
        correction: Procedure,
    ) -> Self {
        Self {
            adj_pvalues: adjust(&pvalues, correction),
            names,
            scores,
            pvalues,
        }
    }

    /// Gets the internal names
    pub fn names(&self) -> &Vec<String> {
        &self.names
    }

    /// Gets the internal scores
    pub fn scores(&self) -> &Vec<f64> {
        &self.scores
    }

    /// Gets the internal pvalues
    pub fn pvalues(&self) -> &Vec<f64> {
        &self.pvalues
    }

    /// Gets the internal adjusted pvalues
    pub fn adj_pvalues(&self) -> &Vec<f64> {
        &self.adj_pvalues
    }

    /// zips the results into a single iterator
    pub fn zip<'a>(&'a self) -> impl Iterator<Item = EntryRRA<'a>> {
        (0..self.names.len()).map(|ix| {
            (
                self.names[ix].as_str(),
                self.scores[ix],
                self.pvalues[ix],
                self.adj_pvalues[ix],
            )
        })
    }
}

#[cfg(test)]
mod testing {

    #[test]
    fn test_results() {
        use super::*;
        use adjustp::Procedure;

        let names = vec!["a".to_string(), "b".to_string()];
        let scores = vec![1.0, 2.0];
        let pvalues = vec![0.1, 0.2];
        let results = ResultsRRA::new(names, scores, pvalues, Procedure::BenjaminiHochberg);

        assert_eq!(results.names().len(), 2);
        assert_eq!(results.scores().len(), 2);
        assert_eq!(results.pvalues().len(), 2);
        assert_eq!(results.adj_pvalues().len(), 2);

        let mut iter = results.zip();
        assert_eq!(iter.next(), Some(("a", 1.0, 0.1, 0.2)));
        assert_eq!(iter.next(), Some(("b", 2.0, 0.2, 0.2)));
        assert_eq!(iter.next(), None);
    }

}

