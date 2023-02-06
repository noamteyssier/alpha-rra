use adjustp::{adjust, Procedure};

type EntryRRA<'a> = (&'a String, f64, f64, f64);

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

    pub fn zip<'a>(&'a self) -> impl Iterator<Item = EntryRRA<'a>> {
        (0..self.names.len()).map(|ix| {
            (
                &self.names[ix],
                self.scores[ix],
                self.pvalues[ix],
                self.adj_pvalues[ix],
            )
        })
    }
}
