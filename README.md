# Opinion Dynamics Incorporating High-Order Interactions

This code accompanies the paper "Opinion Dynamics Incorporating High-Order Interactions". The code is written in *Julia v1.3.1*.

## Datasets

We select a large set of real-world social networks from [Koblenz Network Collection](http://konect.uni-koblenz.de/) and [Network Repository](http://networkrepository.com/). All networks tested are regarded as undirected networks. For those networks that are disconnected originally, we perform our experiments on their largest connected components. 

## Usage

### Network Input

The input of a network consists of the edges in the network. Each line of the input file represents a undirected unweighted edge in the network, which is specified as the format "source_node target_node". Some example files have been put in the `data/` repository.

All input files should be put in the `data/` and have a filename ending with '.txt'. The algorithms will run on these files in lexicographical order of the file names. 

### Run

Execute command `julia main.jl algos hp dist`, where
the directory where the edges list files are.
    
- `algos` denotes which algorithms to run.
    `exact` means running exact algorithm,
    `approx` means running approx algorithm,
    `sparsify` means running the sparsify routine and calcuate the exact solution on the sparsifier,
    and `all` means running all the three algorithms.

- `hp` denotes the choice of hyperparameters.
  `Iter` means testing the algorithm with `[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50, 100, 200, 500, 1000]` iterations,
  `M ` means testing the algorithm with the number of sample edges `M=[1, 10, 100, 200, 500, 1000, 2000]×T×m`,
  `beta` means testing the algorithm with the coefficients equal to `[1 0]` and `[0 1]`, respectively,
  and `best` means testing the algorithm the best hyperparameters. 

- `dist` denotes the distribution of the innate opinions to choose.
  `uni` means uniform distribution,
  `ex` means exponential distribtuion,
  and `pl` means power-law distribtuion.

The result will be printed to both console and file `algos_hp_dist.txt`.

### Examples

```bash
# Clone this repository
$ git clone https://github.com/HODynamic/HODynamic.git

# Go into the source code repository
$ cd src

# Run both the exact and approximate algorithms (ususally for small networks) with power-law distribution
$ julia main.jl all best pl

# Run the approximate algorithms (ususally for large networks) with uniform distribution
$ julia main.jl approx best uni

# Test the sparsification routine with different M
$ julia main.jl sparsify M uni
```
