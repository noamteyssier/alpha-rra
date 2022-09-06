use adjustp::{Procedure, adjust};

pub struct ResultsRRA {
    names: Vec<String>,
    scores: Vec<f64>,
    pvalues: Vec<f64>,
    adj_pvalues: Vec<f64>
}
impl ResultsRRA {
    pub fn new(names: Vec<String>, scores: Vec<f64>, pvalues: Vec<f64>, correction: Procedure) -> Self {
        Self {
            adj_pvalues: adjust(&pvalues, correction),
            names, scores, pvalues,
        }
    }

    pub fn names(&self) -> &Vec<String> {
        &self.names
    }

    pub fn scores(&self) -> &Vec<f64> {
        &self.scores
    }

    pub fn pvalues(&self) -> &Vec<f64> {
        &self.pvalues
    }

    pub fn adj_pvalues(&self) -> &Vec<f64> {
        &self.adj_pvalues
    }
}
