use alpha_rra::alpha_rra;
use ndarray::array;

fn main() {
    let pvalues = array![
        0.3, 0.2, 1e-3, 1e-2, 1e-3, 0.4, 1e-2, 0.4
    ];
    let genes: Vec<String> = vec![
        "gene.0", "gene.0", "gene.1", "gene.1", "gene.1", "gene.2", "gene.2", "gene.2" 
    ].iter().map(|x| x.to_string()).collect();

    let (unique_genes, scores, pvalues) = alpha_rra(&pvalues, &genes, 0.3, 10000);
    println!("{:?}", unique_genes);
    println!("{:?}", scores);
    println!("{:?}", pvalues);
}
