use adjustp::{adjust, Procedure};
use ndarray::Array1;

/// Type alias for the results of an Alpha-RRA run for a single entry
type EntryRRA<'a> = (&'a str, f64, f64, f64);

/// Handles the results of an Alpha-RRA run
pub struct ResultsRRA {
    names: Vec<String>,
    scores: Array1<f64>,
    pvalues: Array1<f64>,
    adj_pvalues: Array1<f64>,
}

impl ResultsRRA {
    /// Creates a new instance
    #[must_use]
    pub fn new(
        names: Vec<String>,
        scores: Array1<f64>,
        pvalues: Array1<f64>,
        correction: Procedure,
    ) -> Self {
        let adj_pvalues = Array1::from_vec(adjust(pvalues.as_slice().unwrap(), correction));
        Self {
            names,
            scores,
            pvalues,
            adj_pvalues,
        }
    }

    /// Gets the internal names
    #[must_use]
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// Gets the internal scores
    #[must_use]
    pub fn scores(&self) -> &Array1<f64> {
        &self.scores
    }

    /// Gets the internal pvalues
    #[must_use]
    pub fn pvalues(&self) -> &Array1<f64> {
        &self.pvalues
    }

    /// Gets the internal adjusted pvalues
    #[must_use]
    pub fn adj_pvalues(&self) -> &Array1<f64> {
        &self.adj_pvalues
    }

    /// zips the results into a single iterator
    pub fn zip(&self) -> impl Iterator<Item = EntryRRA> {
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
    use ndarray::array;

    #[test]
    fn test_results() {
        use super::*;
        use adjustp::Procedure;

        let names = vec!["a".to_string(), "b".to_string()];
        let scores = array![1.0, 2.0];
        let pvalues = array![0.1, 0.2];
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
