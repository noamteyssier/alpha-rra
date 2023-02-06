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
    sizes: &Vec<usize>,
) -> HashMap<usize, Array1<f64>> {
    sizes
        .iter()
        .map(|unique_size| {
            (
                *unique_size,
                run_permutations(n_genes, alpha, npermutations, *unique_size),
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
        let (score, pvalue) = gene_rra(curr, &encode, &nranks, &permutation_vectors, alpha);
        scores[curr] = score;
        pvalues[curr] = pvalue;
    });
    (scores, pvalues)
}

/// Performs the alpha-RRA algorithm
#[must_use]
pub fn alpha_rra(
    pvalues: &Array1<f64>,
    genes: &[String],
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

    // set the number of permutations
    let num_permutations = npermutations * n_genes;

    // calculate rra scores for a vector of random samplings for each unique size
    let permutation_vectors = generate_permutation_vectors(n_genes, alpha, num_permutations, &sizes);

    // calculate empirical pvalues for each of the gene sets given the random nulls
    let (scores, pvalues) =
        calculate_empirical_pvalues(n_genes, &encode, &nranks, &permutation_vectors, alpha);

    // recode the gene names
    let names = recode_index(n_genes, &encode_map);

    // return the results
    ResultsRRA::new(names, scores, pvalues, correction)
}
