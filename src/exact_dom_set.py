"""Radius-constrained clustering algorithms based on minimum dominating
set solvers.
"""

import csv
import os
import shlex
import subprocess
import tempfile
import time
from cProfile import Profile
from pstats import Stats

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

DIR_PATH = os.path.dirname(os.path.realpath(__file__))
EMOS_PATH = os.path.join(
    DIR_PATH,
    "algorithms/EMOS/emos-mds"
)
JAVA_PATH = os.path.abspath(
    os.path.join(
        DIR_PATH,
        "algorithms/MinimumDominatingSet/Code/src/",
    )
)
JAVA_CMD = (
    "java -XX:+ShowCodeDetailsInExceptionMessages "
    f"-cp {JAVA_PATH} heuristic.Main"
)


class MdsClustering:
    """Minimum dominating set-based clustering algorithm.

    Perform clustering using radius constrains using exact or approximate
    minimum dominating set graph algorithms.

    The exact algorithm is provided by Jiang and Zheng. "An exact algorithm
    for the minimum dominating set problem." Proceedings of the Thirty-Second
    International Joint Conference on Artificial Intelligence. 2023.

    The approximate algorithm is provided by Casado et al. "An iterated greedy
    algorithm for finding the minimum dominating set in graphs." Mathematics
    and Computers in Simulation 207 (2023): 41-58.

    Args:
        threshold (float): The threshold T
        manner (str): The manner of clustering between "exact" and "approx"
            Default: "exact"

    Attributes:
        manner (str): The manner of clustering between "exact" and "approx"
        _dist_matrix (numpy.ndarray): The distance matrix
        _adj_matrix (numpy.ndarray): The adjacency matrix
        threshold (float): The threshold T
        clusters (dict): The clusters
        centers (list): The centers of the clusters
        _labels (list): The labels associated with the initial data

    Methods:
        _make_distance_matrix: Compute the distance matrix
        _create_adjacency_matrix: Create the adjacency matrix
        _create_dimacs: Create the dimacs file
        _get_results_exact: Get the results from the exact algorithm
        _get_clusters_exact: Get the clusters from the exact algorithm
        compute_labels: Compute the labels from the clusters
        clustering: Compute the clustering
        plot: Plot the clustering
    """

    def __init__(self, manner="exact"):
        try:
            if not isinstance(manner, str):
                raise TypeError("The manner of clustering must be a string.")
            if manner not in ["exact", "approx"]:
                raise ValueError(
                    "The manner of clustering must be either "
                    "'exact' or 'approx'."
                )
            self._manner = manner
            self._algorithm = (
                self._cluster_exact if manner == "exact" else self._cluster_approx
            )
        except (TypeError, ValueError) as err:
            print(err)
            raise
        self._dist_matrix = None
        self._adj_matrix = None
        self.threshold = None
        self.clusters = None
        self.centers = None
        self._labels = None
        self.exec_time = None
        self._mds_exec_time = None
        self._percentage_of_mds = None
        self.effective_radius = None
        self._temp_dir = None
        self._temp_file = None
        self._result_file = None

    def _make_distance_matrix(self, samples):
        self._dist_matrix = sp.spatial.distance_matrix(samples, samples)

    def _create_adjacency_matrix(self):
        """Converts a distance matrix to an adjacency matrix with edges only for
        distances less than T.
        """

        # Create a boolean mask for the edges to keep
        mask = np.tril(self._dist_matrix <= self.threshold, k=-1)

        # Create the adjacency matrix from the mask
        adjacency_matrix = np.zeros_like(self._dist_matrix)
        adjacency_matrix[mask] = 1
        adjacency_matrix = adjacency_matrix.astype(int)
        self._adj_matrix = adjacency_matrix

    def _create_adj_file(self):
        self._adj_matrix = self._adj_matrix + self._adj_matrix.T
        self._adj_matrix = self._adj_matrix.astype(int)
        with open(self._temp_file, "w") as adj_file:
            adj_file.write(str(self._adj_matrix.shape[0]) + "\n")
            adj_file.write("\n")
            for i in range(self._adj_matrix.shape[0]):
                for j in range(self._adj_matrix.shape[1]):
                    adj_file.write(str(self._adj_matrix[i][j]) + " ")
                adj_file.write("\n")
            adj_file.write("\n")
            adj_file.write("-1")

    def _create_dimacs(self):
        edges = np.sum(self._adj_matrix) // 2
        with open(self._temp_file, "w") as dimacs_file:
            dimacs_file.write(
                "p edge " + str(self._adj_matrix.shape[0]) + " " + str(edges) + "\n"
            )
            for i in range(self._adj_matrix.shape[0]):
                for j in range(i):
                    if self._adj_matrix[i][j] == 1:
                        dimacs_file.write("e " + str(i + 1) + " " + str(j + 1) + "\n")

    def _get_results_exact(self):
        with open(self._result_file, "r") as result_file:
            for line in result_file:
                if line.startswith("Solution:"):
                    self.centers = [
                        int(x) - 1 for x in line.split(":")[1].split(" ") if x.isdigit()
                    ]
                elif line.startswith(">"):
                    self._mds_exec_time = float(line.split(" ")[-1])

    def _get_results_approx(self):
        with open(self._result_file, "r") as result_file:
            reader = csv.reader(result_file, delimiter=";")
            next(reader)
            for row in reader:
                results = row[3].replace("[", "").replace("]", "").split(",")
                self.centers = [int(x) for x in results]
                self._mds_exec_time = float(row[1])

    def _get_clusters_approx(self):
        self.clusters = {}
        for center in self.centers:
            self.clusters[center] = set()
            self.clusters[center].add(center)
        for i in range(self._adj_matrix.shape[0]):
            row = self._dist_matrix[i]
            min_indexes = row[self.centers]
            min_index = np.argmin(min_indexes)
            if self._check_owning(i, self.centers[min_index]):
                self.clusters[self.centers[min_index]].add(i)
            else:
                if "outlier" not in self.clusters:
                    self.clusters["outlier"] = set()
                self.clusters["outlier"].add(i)

    def _check_owning(self, node, center):
        if node == center:
            return True
        if node < center:
            return self._adj_matrix[center][node] == 1
        else:
            return self._adj_matrix[node][center] == 1

    def _get_clusters_exact(self):
        matrix = self._dist_matrix
        dominating_set = self.centers
        clusters = {}
        for node in dominating_set:
            clusters[node] = set()
            clusters[node].add(node)
        for i in range(matrix.shape[0]):
            row = matrix[i]
            min_indexes = row[dominating_set]
            min_index = np.argmin(min_indexes)
            if self._check_owning(i, dominating_set[min_index]):
                clusters[dominating_set[min_index]].add(i)
            else:
                if "outlier" not in clusters:
                    clusters["outlier"] = set()
                clusters["outlier"].add(i)
        self.clusters = clusters

    def _compute_labels(self, samples):
        labels = [0] * (len(samples))
        for i, cluster in enumerate(self.clusters.values()):
            for index in cluster:
                labels[index] = i
        self._labels = labels

    def _cluster_exact(self, samples):
        self._create_dimacs()
        start = time.time()
        result = subprocess.run(
            [EMOS_PATH, self._temp_file],
            capture_output=True,
            text=True,
            check=True,
        )
        end = time.time()
        with open(self._result_file, "w") as result_file:
            result_file.write(result.stdout)
        self._get_results_exact()
        self._get_clusters_exact()
        self._compute_labels(samples)
        self.exec_time = end - start

    def _cluster_approx(self, samples):
        self._create_adj_file()
        cmd = shlex.split(
            JAVA_CMD + f" {self._temp_dir.name} {self._temp_dir.name} {self._temp_file}"
        )
        start = time.time()
        subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
        end = time.time()
        self._get_results_approx()
        self._get_clusters_approx()
        self._compute_labels(samples)
        self.exec_time = end - start

    def _compute_effective_radius(self):
        maximums = []
        for center in self.centers:
            radius = 0
            for i in self.clusters[center]:
                maximums.append(max(radius, self._dist_matrix[center][i]))
        self.effective_radius = max(maximums)

    def _check_dist_matrix(self):
        if self._dist_matrix is None:
            raise ValueError("The distance matrix must be computed before clustering.")

    def _check_adj_matrix(self):
        if self._adj_matrix is None:
            raise ValueError("The adjacency matrix must be computed before clustering.")

    def _init_check(self):
        try:
            self._check_dist_matrix()
            self._check_adj_matrix()
        except ValueError as err:
            print(err)
            raise

    def _check_threshold(self):
        if self.threshold < self.effective_radius:
            raise ValueError("The threshold must be greater than the effective radius.")

    def _check_clusters(self):
        nb_instances = self._dist_matrix.shape[0]
        nb_labels = len(self._labels)
        if nb_instances != nb_labels:
            raise ValueError(
                "The number of labels must be equal to the number of instances."
            )

    def _compute_percentage_of_mds(self):
        self._percentage_of_mds = f"{self._mds_exec_time / self.exec_time * 100}%"

    def _final_check(self):
        try:
            self._check_threshold()
            self._check_clusters()
        except ValueError as err:
            print(err)
            raise

    def clustering(self, samples, threshold):
        """Computer radius-constrained clustering.

        Args:
            samples (NDarray): the population to partition
            threshold (numeric): the radius threshold

        Returns:
            centers: the clusters centers
        """
        self._temp_dir = tempfile.TemporaryDirectory()
        self._temp_file = os.path.join(self._temp_dir.name, "temp.txt")
        self._result_file = os.path.join(self._temp_dir.name, "result.txt")
        self.threshold = threshold
        self._make_distance_matrix(samples)
        self._create_adjacency_matrix()
        self._init_check()
        self._algorithm(samples)
        self._compute_effective_radius()
        self._final_check()
        self._temp_dir.cleanup()
        return self.centers

    def get_exec_time(self):
        return self.exec_time

    def _scores(self):
        outliers_detected = "outlier" in self.clusters
        return (
            self.exec_time,
            self.effective_radius,
            len(self.centers),
            outliers_detected,
        )

    def plot(self, samples):
        fig, ax = plt.subplots()
        sc1 = ax.scatter(samples[:, 0], samples[:, 1], c=self._labels)
        ax.scatter(
            samples[self.centers, 0],
            samples[self.centers, 1],
            c="red",
            label="Centers",
        )
        ax.legend(*sc1.legend_elements(), title="Classes")
        ax.set_title(f"{self._manner} clustering with threshold {self.threshold}")
        plt.show()

    def get_labels(self, samples):
        if self._labels is None:
            self._compute_labels(samples)
        return self._labels

    def __str__(self) -> str:
        return f"MdsClustering(threshold={self.threshold}, manner={self._manner})"

    def __repr__(self) -> str:
        return f"MdsClustering(threshold={self.threshold}, manner={self._manner})"

    def _profile(self, samples):
        profiler = Profile()
        profiler.runcall(self.clustering, samples)
        stats = Stats(profiler)
        stats.strip_dirs()
        stats.sort_stats("cumulative")
        stats.print_stats()
