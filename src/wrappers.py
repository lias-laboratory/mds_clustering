"""Python wrappers for equiwide and protoclust."""

import os
import subprocess
import tempfile
import time

import numpy as np
from scipy.spatial import distance_matrix

from algorithms.ewclustering.equiwide_clustering import EquiwideClustering

os.environ["R_LIBS_USER"] = os.path.join(os.path.dirname(__file__), "algorithms")
PROTOCLUST = os.path.join(os.path.dirname(__file__), "algorithms/protoclust.R")


class EQWClustering:
    """Python wrapper for equiwide."""

    def __init__(self, cover="lp"):
        self._eqw = EquiwideClustering(
            measure="radius",
            cover=cover,
            verbose=0,
        )
        self.exec_time = None

    def compute_centers(self, dist, labels):
        """Compute centers from labels."""
        centers = []
        for label in np.unique(labels):
            points = np.where(labels == label)[0]
            max_indexes = []
            for point in points:
                max_dist_index = np.argmax(dist[point][points])
                max_indexes.append(points[max_dist_index])
            min_index = np.argmin(dist[max_indexes, points])
            centers.append(points[min_index])
        return centers

    def clustering(self, samples, threshold):
        """Computer radius-constrained clustering based on equiwide."""
        dist = distance_matrix(samples, samples)
        start_time = time.time()
        labels = self._eqw.clustering(dist, local_threshold=threshold)
        self.exec_time = time.time() - start_time
        return self.compute_centers(dist, labels)

    def get_exec_time(self):
        """Retrieve core algorithm execution time."""
        return self.exec_time


class Protoclust:
    """Python wrapper for protoclust."""

    def __init__(self):
        self.exec_time = None

    def clustering(self, samples, threshold):
        """Computer radius-constrained clustering based on protoclust."""
        dist = distance_matrix(samples, samples)
        with tempfile.NamedTemporaryFile() as tmp:
            np.savetxt(tmp.name, dist, delimiter=",")
            results = subprocess.run(
                [PROTOCLUST, tmp.name, str(threshold)],
                capture_output=True,
                text=True,
                check=True,
            )
        duration, *centers = results.stdout.split(",")
        self.exec_time = float(duration)
        centers = [int(center) for center in centers]
        return centers

    def get_exec_time(self):
        """Retrieve core algorithm execution time."""
        return self.exec_time
