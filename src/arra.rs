use crate::{utils::recode_index, ResultsRRA};
use adjustp::Procedure;
use anyhow::{bail, Result};
use getset::Getters;
use hashbrown::HashMap;
use ndarray::Array1;

use super::{
    filter_alpha, group_sizes, normed_ranks,
    permutations::run_permutations,
    robust_rank::robust_rank_aggregation,
    utils::{empirical_cdf, encode_index, select_ranks},
};

/// Calculates an empirical p-value of the robust rank aggregation for the current gene set with
/// respect to random permutations of that size
///
/// # Arguments
/// * `current_idx` - The current index (i.e. gene index)
/// * `encodings` - The encodings (i.e. gene indices)
/// * `nranks` - The normalized ranks
/// * `permutation_vectors` - The permutation vectors
/// * `alpha` - The alpha threshold value
///
/// # Returns
/// A tuple of the score and the empirical p-value
fn gene_rra(
    current_idx: usize,
    encodings: &[usize],
    nranks: &Array1<f64>,
    permutation_vectors: &HashMap<usize, Array1<f64>>,
    alpha: f64,
) -> (f64, f64) {
    let gene_ranks = select_ranks(current_idx, encodings, nranks);
    let filtered = filter_alpha(&gene_ranks, alpha);
    let score = robust_rank_aggregation(&filtered, gene_ranks.len());
    (
        score,
        empirical_cdf(
            score,
            permutation_vectors
                .get(&gene_ranks.len())
                .expect("Unexpected missing key"),
        ),
    )
}

/// Generates a vector of random samplings for each unique size then scores them
///
/// # Arguments
/// * `n_genes` - The number of genes
/// * `alpha` - The alpha threshold value
/// * `npermutations` - The number of permutations
/// * `sizes` - The unique sizes of each group
///
/// # Returns
/// A map of the unique sizes to the random samplings scores
fn generate_permutation_vectors(
    n_genes: usize,
    alpha: f64,
    npermutations: usize,
    sizes: &[usize],
    seed: u64,
) -> HashMap<usize, Array1<f64>> {
    sizes
        .iter()
        .map(|unique_size| {
            (
                *unique_size,
                run_permutations(n_genes, alpha, npermutations, *unique_size, seed),
            )
        })
        .map(|(u, v)| (u, Array1::from_vec(v)))
        .collect::<HashMap<usize, Array1<f64>>>()
}

/// Calculates the empirical p-values for each of the gene sets given the random nulls
///
/// # Arguments
/// * `n_genes` - The number of genes
/// * `encode` - The encodings (i.e. gene indices)
/// * `nranks` - The normalized ranks
/// * `permutation_vectors` - The permutation vectors
/// * `alpha` - The alpha threshold value
///
/// # Returns
/// A tuple of the scores and the empirical p-values
fn calculate_empirical_pvalues(
    n_genes: usize,
    encode: &[usize],
    nranks: &Array1<f64>,
    permutation_vectors: &HashMap<usize, Array1<f64>>,
    alpha: f64,
) -> (Array1<f64>, Array1<f64>) {
    let mut scores = Array1::zeros(n_genes);
    let mut pvalues = Array1::zeros(n_genes);
    (0..n_genes).for_each(|curr| {
        let (score, pvalue) = gene_rra(curr, encode, nranks, permutation_vectors, alpha);
        scores[curr] = score;
        pvalues[curr] = pvalue;
    });
    (scores, pvalues)
}

/// A struct to perform the Alpha RRA Algorithm
#[derive(Getters)]
pub struct AlphaRRA {
    /// A map of the gene names to their index
    #[getset(get = "pub")]
    encode_map: HashMap<usize, String>,

    /// The encodings (i.e. gene indices)
    #[getset(get = "pub")]
    encode: Vec<usize>,

    /// The alpha threshold value
    #[getset(get_copy = "pub")]
    alpha: f64,

    /// The number of permutations
    #[allow(dead_code)]
    #[getset(get_copy = "pub")]
    n_permutations: usize,

    /// The correction method
    #[getset(get_copy = "pub")]
    correction: Procedure,

    /// The number of unique genes
    #[getset(get_copy = "pub")]
    n_genes: usize,

    /// The permutation vectors (i.e. random samplings and their alpha RRA scores)
    #[getset(get = "pub")]
    permutation_vectors: HashMap<usize, Array1<f64>>,

    /// The random seed to use in generating permutation vectors
    #[allow(dead_code)]
    #[getset(get_copy = "pub")]
    seed: u64,
}
impl AlphaRRA {
    /// Creates a new AlphaRRA struct
    ///
    /// # Arguments
    /// * `genes` - The gene names
    /// * `alpha` - The alpha threshold value
    /// * `n_permutations` - The number of permutations
    /// * `correction` - The correction method
    pub fn new(
        genes: &[String],
        alpha: f64,
        n_permutations: usize,
        correction: Procedure,
        seed: u64,
    ) -> Self {
        let (encode_map, encode) = encode_index(genes);
        let n_genes = encode_map.len();
        let n_permutations = n_permutations * n_genes;
        let permutation_vectors = generate_permutation_vectors(
            n_genes,
            alpha,
            n_permutations,
            &group_sizes(&encode),
            seed,
        );
        Self {
            encode_map,
            encode,
            alpha,
            n_permutations,
            correction,
            n_genes,
            permutation_vectors,
            seed,
        }
    }

    /// Runs the Alpha RRA Algorithm for a given set of p-values
    ///
    /// # Arguments
    /// * `pvalues` - The p-values
    ///
    /// # Returns
    /// A ResultsRRA struct
    pub fn run(&self, pvalues: &Array1<f64>) -> Result<ResultsRRA> {
        if pvalues.len() != self.encode.len() {
            bail!("The number of p-values does not match the number of sgRNAs in the library");
        }
        let nranks = normed_ranks(pvalues);
        let (scores, pvalues) = calculate_empirical_pvalues(
            self.n_genes,
            &self.encode,
            &nranks,
            &self.permutation_vectors,
            self.alpha,
        );
        let names = recode_index(self.n_genes, &self.encode_map);
        let result = ResultsRRA::new(names, scores, pvalues, self.correction);
        Ok(result)
    }
}
