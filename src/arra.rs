use crate::{utils::recode_index, ResultsRRA};
use adjustp::Procedure;
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

/// Performs the alpha-RRA algorithm
#[must_use]
pub fn alpha_rra(
    pvalues: &Array1<f64>,
    genes: &Vec<String>,
    alpha: f64,
    npermutations: usize,
    correction: Procedure,
) -> ResultsRRA {
    // encode the gene names
    let (encode_map, encode) = encode_index(genes);
    let n_genes = encode_map.len();

    // calculate the normalized ranks
    let nranks = normed_ranks(pvalues);

    // group the sizes of the gene sets
    let sizes = group_sizes(&encode);

    // calculate rra scores for a vector of random samplings for each unique size
    let permutation_vectors = sizes
        .iter()
        .map(|unique_size| {
            (
                *unique_size,
                run_permutations(nranks.len(), alpha, npermutations * n_genes, *unique_size),
            )
        })
        .map(|(u, v)| (u, Array1::from_vec(v)))
        .collect::<HashMap<usize, Array1<f64>>>();

    // calculate empirical pvalues for each of the gene sets given the random nulls
    let (scores, pvalues): (Vec<f64>, Vec<f64>) = (0..n_genes)
        .map(|curr| gene_rra(curr, &encode, &nranks, &permutation_vectors, alpha))
        .unzip();

    // recode the gene names
    let names = recode_index(n_genes, &encode_map);

    // return the results
    ResultsRRA::new(names, scores, pvalues, correction)
}
