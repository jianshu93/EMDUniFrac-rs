use anyhow::{Context, Result};
use clap::{Arg, Command};
use phylotree::tree::Tree;
use rayon::prelude::*;
use std::{
    collections::HashSet,
    fs::File,
    io::{BufRead, BufReader, Write},
    path::Path,
};

static mut MASTER_L_TOTAL: f64 = 0.0;

fn main() -> Result<()> {
    // Initialize logger
    println!("\n ************** initializing logger *****************\n");
    let _ = env_logger::Builder::from_default_env().init();
    let matches = Command::new("Unweighted_UniFrac")
        .version("0.1.0")
        .about("Fast Unweighted UniFrac using EMDUniFrac")
        .arg(Arg::new("tree")
            .short('t')
            .long("tree")
            .value_name("TREE_FILE")
            .help("Input newick format tree file")
            .required(true))
        .arg(Arg::new("table")
            .short('i')
            .long("input")
            .value_name("TABLE_FILE")
            .help("Input tab-delimited sample-feature table")
            .required(true))
        .arg(Arg::new("output")
            .short('o')
            .long("output")
            .value_name("OUTPUT_FILE")
            .help("Output file for distance matrix")
            .required(true))
        .get_matches();

    let tree_file = matches.get_one::<String>("tree").unwrap();
    let table_file = matches.get_one::<String>("table").unwrap();
    let output_file = matches.get_one::<String>("output").unwrap();

    // Read the tree
    let tree = Tree::from_file(Path::new(tree_file))?;

    // Read and normalize the sample-feature table
    let (taxa_order, sample_names, mut presence_matrix) = read_sample_table(table_file)?;
    normalize_samples(&mut presence_matrix);
    let n_samples = sample_names.len();

    // Build Tint and lint from the given tree (no pruning)
    let (tint, lint, nodes_in_order, node_name_map) = build_tint_lint(&tree)?;

    // Compute MASTER_L_TOTAL: sum of all branch lengths
    let master_l_total: f64 = lint.iter().sum();
    unsafe {
        MASTER_L_TOTAL = master_l_total;
    }

    // Map leaf names to node indices
    let leaf_map = build_leaf_map(&tree, &node_name_map)?;

    // Prepare all pairs i,j
    let pairs: Vec<(usize, usize)> = (0..n_samples)
        .flat_map(|i| (i+1..n_samples).map(move |j| (i, j)))
        .collect();

    // Compute all pairs in parallel
    let updates: Vec<(usize, usize, f64)> = pairs.par_iter().map(|&(i, j)| {
        let uni = compute_unifrac_for_pair(
            &tint,
            &lint,
            &nodes_in_order,
            &leaf_map,
            &taxa_order,
            &presence_matrix,
            i,
            j
        ).unwrap();
        (i, j, uni)
    }).collect();

    // Construct distance matrix
    let mut dist_matrix = vec![0.0; n_samples * n_samples];
    for i in 0..n_samples {
        dist_matrix[i * n_samples + i] = 0.0;
    }

    for &(i, j, uni) in &updates {
        dist_matrix[i * n_samples + j] = uni;
        dist_matrix[j * n_samples + i] = uni;
    }

    // Write output matrix
    write_matrix(&sample_names, &dist_matrix, n_samples, output_file)?;

    Ok(())
}

/// Read the sample-feature table.
/// Any value > 0 is converted to 1.0, else 0.0 for unweighted presence/absence.
fn read_sample_table(filename: &str) -> Result<(Vec<String>, Vec<String>, Vec<Vec<f64>>)> {
    let f = File::open(filename)?;
    let mut lines = BufReader::new(f).lines();

    // First line: parse sample names
    let header = lines.next().context("No header in table")??;
    let mut hdr_split = header.split_whitespace();
    hdr_split.next(); // ignore the first element
    let sample_names: Vec<String> = hdr_split.map(|s| s.to_string()).collect();

    let mut taxa_order = Vec::new();
    let mut presence_matrix = Vec::new();

    for line in lines {
        let line = line?;
        let mut parts = line.split_whitespace();
        let taxon = parts.next().context("Taxon missing in a line")?.to_string();
        taxa_order.push(taxon);
        let values: Vec<f64> = parts.map(|x| {
            let val: f64 = x.parse().unwrap_or(0.0);
            if val > 0.0 {1.0} else {0.0}
        }).collect();
        presence_matrix.push(values);
    }

    Ok((taxa_order, sample_names, presence_matrix))
}

/// Normalize each sample so that they form a probability distribution (sum to 1).
fn normalize_samples(presence_matrix: &mut [Vec<f64>]) {
    let n_samples = presence_matrix[0].len();
    for s in 0..n_samples {
        let sum: f64 = presence_matrix.iter().map(|row| row[s]).sum();
        if sum > 0.0 {
            for row in presence_matrix.iter_mut() {
                row[s] /= sum;
            }
        }
    }
}

/// Build Tint and lint from the given tree:
fn build_tint_lint(tree: &Tree) -> Result<(Vec<usize>, Vec<f64>, Vec<usize>, Vec<Option<String>>)> {
    let root = tree.get_root()?;
    let postord = tree.postorder(&root)?;
    let num_nodes = postord.len();

    let mut pos_map = vec![0; tree.size()];
    for (i, &nid) in postord.iter().enumerate() {
        pos_map[nid] = i;
    }

    let mut tint = vec![0; num_nodes];
    let mut lint = vec![0.0; num_nodes];
    let mut node_name_map = vec![None; num_nodes];

    let root_idx = pos_map[root];
    tint[root_idx] = root_idx;
    lint[root_idx] = 0.0;

    for &nid in &postord {
        let node = tree.get(&nid)?;
        if let Some(name) = &node.name {
            node_name_map[pos_map[nid]] = Some(name.clone());
        }
        if nid != root {
            let p = node.parent.ok_or_else(|| anyhow::anyhow!("Node has no parent but is not root"))?;
            tint[pos_map[nid]] = pos_map[p];
            lint[pos_map[nid]] = node.parent_edge.unwrap_or(0.0);
        }
    }

    Ok((tint, lint, postord, node_name_map))
}

/// Build a map from leaf name to its index in nodes_in_order
fn build_leaf_map(
    tree: &Tree,
    node_name_map: &[Option<String>]
) -> Result<std::collections::HashMap<String, usize>> {
    let leaves = tree.get_leaves();
    let mut leaf_map = std::collections::HashMap::new();

    let root = tree.get_root()?;
    let postord = tree.postorder(&root)?;
    let mut pos_map = vec![0; tree.size()];
    for (i, &nid) in postord.iter().enumerate() {
        pos_map[nid] = i;
    }

    for l in leaves {
        let node = tree.get(&l)?;
        if node.is_tip() {
            if let Some(lname) = &node.name {
                let idx = pos_map[l];
                leaf_map.insert(lname.clone(), idx);
            }
        }
    }

    Ok(leaf_map)
}

/// Compute UniFrac for a given pair of samples i,j using the EMDUnifrac_unweighted logic:
/// - partial_sums at leaves = P - Q
/// - Propagate up the tree, accumulate Z
/// - UniFrac = Z
fn compute_unifrac_for_pair(
    tint: &[usize],
    lint: &[f64],
    nodes_in_order: &[usize],
    leaf_map: &std::collections::HashMap<String, usize>,
    taxa_order: &[String],
    presence_matrix: &[Vec<f64>],
    i: usize,
    j: usize
) -> Result<f64> {
    let num_nodes = nodes_in_order.len();
    let mut partial_sums = vec![0.0; num_nodes];

    // Compute leaf differences in parallel
    let diffs = taxa_order.par_iter()
        .enumerate()
        .filter_map(|(t_idx, taxon)| {
            if let Some(&leaf_idx) = leaf_map.get(taxon) {
                let val_i = presence_matrix[t_idx][i];
                let val_j = presence_matrix[t_idx][j];
                // If no difference, skip
                if (val_i - val_j).abs() > 1e-14 {
                    Some((leaf_idx, val_i - val_j))
                } else {
                    None
                }
            } else {
                None
            }
        })
        .collect::<Vec<(usize,f64)>>();

    for (leaf_idx, diff) in diffs {
        partial_sums[leaf_idx] = diff;
    }

    // Propagate differences up the tree
    let mut Z = 0.0;
    for node_pos in 0..num_nodes {
        if tint[node_pos] == node_pos {
            continue; // root
        }
        let val = partial_sums[node_pos];
        partial_sums[tint[node_pos]] += val;
        Z += lint[node_pos] * val.abs();
    }
    let unifrac = Z;
    Ok(unifrac)
}

/// Write the resulting matrix to a file
fn write_matrix(sample_names: &[String], dist_matrix: &[f64], n: usize, output_file: &str) -> Result<()> {
    let mut file = File::create(output_file)?;
    write!(file, "Sample")?;
    for sn in sample_names {
        write!(file, "\t{}", sn)?;
    }
    writeln!(file)?;

    for i in 0..n {
        write!(file, "{}", sample_names[i])?;
        for j in 0..n {
            write!(file, "\t{:.6}", dist_matrix[i * n + j])?;
        }
        writeln!(file)?;
    }

    Ok(())
}
