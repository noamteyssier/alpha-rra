use super::{filter_alpha, robust_rank_aggregation};
use ndarray::Array1;
use ndarray_rand::{
    rand::{Rng, SeedableRng},
    rand_distr::Uniform,
};
use rand_chacha::ChaCha8Rng;
use rayon::prelude::*;

/// Sample a vector of normalized ranks given the number of ranks and the number of samples
///
/// # Arguments
/// * `num_ranks` - The number of ranks
/// * `num_samples` - The number of samples
///
/// # Returns
/// A vector of normalized ranks
fn sample_normed_ranks(num_ranks: usize, num_samples: usize, idx: usize, seed: u64) -> Array1<f64> {
    let rng = ChaCha8Rng::seed_from_u64(seed + idx as u64);
    rng.sample_iter(Uniform::new(1usize, num_ranks + 1))
        .take(num_samples)
        .map(|x| x as f64 / num_ranks as f64)
        .collect()
}

/// Perform the robust rank aggregation algorithm on a permuted set of normalized ranks
///
/// # Arguments
/// * `num_ranks` - The number of ranks
/// * `alpha` - The alpha value
/// * `npermutations` - The number of permutations
/// * `unique_size` - The number of unique values
///
/// # Returns
/// A vector of the robust rank aggregation values
#[must_use]
pub fn run_permutations(
    num_ranks: usize,
    alpha: f64,
    npermutations: usize,
    unique_size: usize,
    seed: u64,
) -> Vec<f64> {
    (0..npermutations)
        .into_par_iter()
        .map(|idx| sample_normed_ranks(num_ranks, unique_size, idx, seed))
        .map(|choices| filter_alpha(&choices, alpha))
        .map(|filtered| robust_rank_aggregation(&filtered, unique_size))
        .collect()
}

#[cfg(test)]
mod testing {
    use super::run_permutations;
    use ndarray::Array1;
    use ndarray_rand::rand_distr::Uniform;
    use ndarray_rand::RandomExt;
    use statrs::statistics::Statistics;

    #[test]
    fn test_run() {
        let unique_size = 5;
        let num_samples = 100;
        let alpha = 0.3;
        let npermutations = 1000;
        let permutations = run_permutations(num_samples, alpha, npermutations, unique_size, 0);
        assert_eq!(permutations.len(), npermutations);
    }

    #[test]
    fn test_sample_int() {
        for _ in 0..100 {
            let num_samples = 100;
            let size = 5;
            let nranks = Array1::random((size,), Uniform::new(1, num_samples))
                .iter()
                .map(|x| f64::from(*x) / f64::from(num_samples))
                .collect::<Array1<f64>>();
            assert_eq!(nranks.len(), size);
            assert!((&nranks).max() <= 1.);
            assert!(nranks.min() > 0.);
        }
    }

    #[test]
    fn seeded_normed_ranks() {
        let num_samples = 100;
        let size = 5;
        let seed = 0;
        let ranks_1 = super::sample_normed_ranks(num_samples, size, 0, seed);
        for _ in 0..100 {
            let ranks_2 = super::sample_normed_ranks(num_samples, size, 0, seed);
            assert_eq!(ranks_1, ranks_2);
        }
        for _ in 0..100 {
            let ranks_3 = super::sample_normed_ranks(num_samples, size, 1, seed);
            assert_ne!(ranks_1, ranks_3);
        }
    }
}
