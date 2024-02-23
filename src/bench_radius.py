"""A clustering under radius constraint solutions benchmark.

This script runs four clustering algorithms:
  - EMOS (https://github.com/huajiang-ynu/ijcai23-mds/)
  - A MinimumDominatingSet heuristic (https://github.com/AlejandraCasado/MinimumDominatingSet/)
  - Protoclust (https://cran.r-project.org/web/packages/protoclust/)
  - Equiwide-clustering

Experimental results are stored in a csv file (results/bench_results.csv).

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
"""

import os
import sys
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from pathlib import Path

import pandas as pd
from scipy.spatial import distance_matrix
from sklearn.datasets import fetch_openml
from tqdm import tqdm

import requirements
from exact_dom_set import MdsClustering
from wrappers import EQWClustering, Protoclust

DATASETS_FAST = [
    ("iris", 1.4282856857085696),
    ("wine", 232.083),
    ("glass", 3.94),
    ("ionosphere", 5.4524),
    ("wdbc", 1197.42),
    ("synthetic_control", 70.112),
    ("vehicle", 155.046),
]
DATASETS_SLOW = [
    ("yeast", 0.42296571965113201),
    ("ozone-level-8hr", 245.5827),
    ("waveform-5000", 10.732427498008081),
]
DATASETS = {
    "fast": DATASETS_FAST,
    "slow": DATASETS_SLOW,
    "all": DATASETS_FAST + DATASETS_SLOW,
}

ROOT_PATH = os.path.dirname(Path(os.path.realpath(__file__)).parent)
RAW_CSV = os.path.join(ROOT_PATH, "results/raw_results.csv")
BENCH_CSV = os.path.join(ROOT_PATH, "results/bench_results.csv")


def export_results(results):
    """Export raw and aggregated results to csv.

    Args:
        results (DataFrame): the raw results DataFrame
    """
    results.to_csv(RAW_CSV, index=False)
    ordering = results["dataset"].drop_duplicates()
    (
        results.groupby(["dataset", "algorithm"])
        .agg({
            "effective_radius": ["mean", "std"],
            "nb_classes": ["min", "max"],
            "exec_time": ["mean", "std"],
        })
        .unstack()
        .swaplevel(i=0, j=2, axis=1)
        .swaplevel(i=1, j=2, axis=1)
        .sort_index(axis=1)
        .loc[ordering]
        .reset_index()
        .to_csv(BENCH_CSV, index=False)
    )


def compute_effective_radius(dist, centers):
    """Compute the maximum radius.

    Args:
        dist: the distance matrix
        centers (list): the centers indexes

    Returns:
        max_radius (float): the maximum radius
    """
    return dist[:, centers].min(axis=1).max()


def compute(datasets, nb_repetitions, disable_progress=False):
    """Run clustering algorithms and store their results.

    Args:
        datasets (list): list of dataset names and thresholds to use for the benchmark
        nb_repetitions (int): number of repetitions

    Returns:
        results (DataFrame): the raw results DataFrame
    """
    results = pd.DataFrame()
    for dataset, threshold in tqdm(
        datasets,
        position=0,
        desc="datasets",
        disable=disable_progress,
    ):
        data = fetch_openml(name=dataset, version=1, parser="auto").data
        dist = distance_matrix(data, data)
        for cls, algorithm in [
            (EQWClustering(), "equiwide"),
            (Protoclust(), "protoclust"),
            (MdsClustering(manner="exact"), "mds-exact"),
            (MdsClustering(manner="approx"), "mds-approx"),
        ]:
            for i in tqdm(
                range(nb_repetitions),
                position=1,
                desc=f"{dataset} ({algorithm})",
                leave=False,
                disable=disable_progress,
            ):
                centers = cls.clustering(data, threshold)
                score = pd.DataFrame(
                    data=[{
                        "dataset": dataset,
                        "algorithm": algorithm,
                        "iteration": i + 1,
                        "effective_radius": compute_effective_radius(dist, centers),
                        "nb_classes": len(centers),
                        "exec_time": cls.get_exec_time(),
                    }]
                )
                results = pd.concat([results, score])
    return results.reset_index(drop=True)


def check_algorithms():
    """Run a quick benchmark."""
    compute(DATASETS_FAST[0:1], 1, disable_progress=True)


def main():
    """Parse command line arguments and run the benchmark."""
    parser = ArgumentParser(
        description=__doc__,
        epilog="Example: python3 bench_radius.py -d fast -n 5",
        formatter_class=RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-d",
        "--dataset",
        choices=["fast", "slow", "all"],
        default="fast",
        help="The subset of datasets to use (fast, slow or all). Default: fast",
    )
    parser.add_argument(
        "-n",
        "--nb-repetitions",
        type=int,
        default=10,
        help="The number of repetitions of the clustering algorithm. Default: 10",
    )
    parser.add_argument(
        "-c",
        "--check",
        action="store_true",
        help="Check if the installation is correct and tools are installed",
    )
    parser.add_argument(
        "--check-requirements",
        action="store_true",
        help="Check if the requirements are fulfilled",
    )
    args = parser.parse_args()
    if args.check:
        requirements.check_requirements()
        check_algorithms()
        sys.exit(0)
    if args.check_requirements:
        requirements.check_requirements()
        sys.exit(0)
    scores = compute(DATASETS[args.dataset], args.nb_repetitions)
    export_results(scores)


if __name__ == "__main__":
    main()
