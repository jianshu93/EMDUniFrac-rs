# EMDUniFrac implementation in Rust

This is an example repo to show how to compute the Earth Mover Distace-based UniFrac distance in Rust. 
It uses the [phylotree](https://github.com/lucblassel/phylotree-rs) crate to parse the tree file and feature-sample table (OTU table for example) and then compute the EMDUniFrac distance. The code is full parallelized for many samples.




## Install
Install Rustup (Rust toolchain) here first: https://rustup.rs
```bash
git clone https://github.com/jianshu93/EMDUniFrac-rs
cd EMDUniFrac-rs
cargo build --release
./target/release/EMDUniFrac -h
```

## Usage 
```bash
 ************** initializing logger *****************

Fast Earth Mover Distance (EMD) UniFrac

Usage: EMDUniFrac [OPTIONS] --tree <TREE_FILE> --input <TABLE_FILE> --output <OUTPUT_FILE>

Options:
  -t, --tree <TREE_FILE>      Input newick format tree file
  -i, --input <TABLE_FILE>    Input tab-delimited sample-feature table
  -o, --output <OUTPUT_FILE>  Output file for distance matrix
      --weighted              Weighted EMDUniFrac
  -h, --help                  Print help
  -V, --version               Print version
```

### example
```bash
### remove bootstrap support first if you have it

### Then run unifrac like this:
./target/release/EMDUniFrac -t data/test_rot_new1.nwk -i data/table.txt -o try.txt
cat try.txt
```

## References
1.McClelland, J. and Koslicki, D., 2018. EMDUniFrac: exact linear time computation of the UniFrac metric and identification of differentially abundant organisms. Journal of mathematical biology, 77, pp.935-949.
