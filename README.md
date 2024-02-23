# Feasibility report dominating set

This repository contains the source code and reproducibility instructions for the feasibility study about minimum dominating set implementation of clustering under radius constraints.

## Installation

### Requirements

- Python 3 (tested with 3.11)
- Java Development Kit (tested with OpenJDK 17)
- R (tested with 4.2.2)
- GCC (tested with 12.2.0)
- Git

### Getting the code

The source code is available on [github](https://github.com/lias-laboratory/mds_clustering):

```bash
git clone https://github.com/lias-laboratory/mds_clustering
cd mds_clustering
```

### Quick install

Once you have verified that you have the required software (see the requirements section above), you can install the additional dependencies by using the provided setup script:

```bash
./setup.sh
```

If you prefer to install and configure dependencies manually, please refer to the INSTALL.md file.

## Usage

### Benchmarking

The `bench_radius.py` file is the entry point for the benchmark, from which you can run experiments.
The following instruction runs the benchmark using the article settings (all datasets, 10 iterations).

```bash
python3 src/bench_radius.py -d all
```

If you need help with the command line arguments, please refer to the help option of the benchmark.

```bash
python3 src/bench_radius.py --help
```

## Comparing results

Upon completion, benchmark results are saved in two files:

- `results/raw_results.csv`, which contains results for each dataset, algorithm and iteration,
- `results/bench_results.csv`, which contains results aggregated by dataset and algorithm.

## License

This project in licensed under the MIT License. You can find the full text of the license in the `LICENSE` file.

## Contributors

This program is part of the research work on the clustering under radius constraint problem, with the support of the LIAS laboratory.

<p align="left">
<a href=https://www.lias-lab.fr/fr/><img src=https://www.lias-lab.fr/images/logo_lias.png width=200></a>
</p>

This project was made by the following contributors:

- [Quentin Haenn](https://www.lias-lab.fr/members/quentinhaenn/) (Main contributor)
- [Brice Chardin](https://www.lias-lab.fr/members/bricechardin/) (Reviewer and supervisor)
- [Mickael Baron](https://www.lias-lab.fr/members/mickaelbaron/) (Reviewer and supervisor)

## Acknowledgments

We would like to thank the following people for kindly providing us with the source code of their algorithms:

- [Alejandra Casado](https://github.com/AlejandraCasado) for the minimum dominating set heuristic code, presented in [An iterated greedy algorithm for finding the minimum dominating set in graphs](https://www.sciencedirect.com/science/article/pii/S0378475422005055)
- [Hua Jiang](https://github.com/huajiang-ynu) for the minimum dominating set exact algorithm code presented in [An exact algorithm for the minimum dominating set problem](https://dl.acm.org/doi/abs/10.24963/ijcai.2023/622)
