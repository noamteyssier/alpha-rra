pub mod arra;
pub mod permutations;
pub mod results;
pub mod robust_rank;
pub mod utils;

pub use arra::AlphaRRA;
pub use results::ResultsRRA;
use robust_rank::robust_rank_aggregation;
use utils::{filter_alpha, group_sizes, normed_ranks};
