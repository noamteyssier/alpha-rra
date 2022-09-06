pub mod arra;
pub mod robust_rank;
pub mod utils;
pub mod permutations;
pub mod results;

pub use arra::alpha_rra;
pub use results::ResultsRRA;
use robust_rank::robust_rank_aggregation;
use utils::{normed_ranks, group_sizes, filter_alpha};
