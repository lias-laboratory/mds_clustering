#!/usr/bin/env python3
"""Equiwide-clustering : Clustering algorithm based on a maximum intra-cluster distance.

This module entry point is the clustering class with the clustering method.
"""

import numpy as np
from mip import BINARY, Model, OptimizationStatus, minimize, xsum


class EquiwideClustering:
    """Class for a clustering algorithm based on a maximum intra-cluster
    distance. Generate clusters from a distance matrix.

    Attributes:
        _measure (string): Which distance is considered between radius and
            diameter
        _cover (string): Which method is used to compute minimal set cover
        _threshold (float): The distance threshold, i.e. the maximum diameter
            of clusters for the max-diameter algorithm or the maximum distance
            to the cluster representative for the max-radius algorithm
        _cutoff (float): When cover is lp, a cutoff can be set to limit runtime
            of the computation of minimal set cover
        _solver (str): The name of the solver to use with MIP
        _centers (list): The list of the centers associated with each
            compatible sets on radius
        _verbose (int): Verbose : control which messages are printed
    """

    def __init__(
        self,
        measure="radius",
        cover="greedy",
        threshold=None,
        cutoff=None,
        solver=None,
        verbose=1,
    ):
        """Constructor

        Args:
            measure (str): {"radius"} to choose which algorithm use. Currently
                only "radius" is supported
            cover (str): {"greedy", "lp"} to choose how is computed the min
                set cover. "greedy" for an approximate greedy algorithm. "lp"
                for a linear programming algorithm. Defaults to "greedy"
            threshold (float): The distance threshold, i.e. the maximum
                diameter of clusters for the max-diameter algorithm or the
                maximum distance to the cluster representative for the
                max-radius algorithm. Threshold can be set later. Defaults
                to None
            cutoff (float): When cover is lp, a cutoff can be set to limit
                runtime of the computation of minimal set cover. cutoff in
                seconds. Defaults to None
            solver (str): The name of the solver to use with MIP: typically GRB
                or CBC. Defaults to None: the solver is automatically chosen:
                Gurobi if available, else CBC
            verbose (int): 0 to disable all messages printed, 1 to enable
                messages about current clustering, 2 to enable messages to
                benchmark the algorithm

        Raises:
            ValueError: When measure parameter is not "diameter" or "radius".
                Or when the cover parameter is not "greedy" or "lp". Or when a
                cutoff is set with another cover than "lp".
            TypeError: If threshold or cutoff is not None and is not and
                integer or a float.
        """
        # measure should be either radius or diameter
        if measure in {"diameter", "radius"}:
            self._measure = measure
        else:
            raise ValueError("measure should be either 'radius' or 'diameter'")
        # cover should be either greedy, either lp, or exact
        if cover in {"greedy", "lp", "exact"}:
            self._cover = cover
        else:
            raise ValueError("cover should be 'greedy', 'exact' or 'lp'")
        # threshold should be an integer or a float
        if threshold is None:
            self._threshold = None
        else:
            if type(threshold) in {int, float}:
                self._threshold = threshold
            else:
                raise TypeError("threshold should be a number (integer or float)")
        # cutoff should be an integer or a float
        if cutoff is None:
            self._cutoff = None
        else:
            if self._cover != "lp":
                raise ValueError("cutoff can only be set with 'lp' for cover")
            if type(cutoff) in {int, float}:
                self._cutoff = cutoff
            else:
                raise TypeError("cutoff should be a number (integer or float)")
        self._solver = solver
        if verbose == 1:
            self._verbose = 1
        elif verbose > 1:
            self._verbose = 2
        else:
            self._verbose = 0
        self._centers = None

    def _compute_distances(self, samples, dist_function):
        """Compute the distance matrix from samples data using the distance
        function provided.

        Compute the distance matrix from the lines of samples and
        using the dist_function to compute the distance between two lines.

        Args:
            samples (np.darray or list): An array of samples where each ligne
                is a sample and each column is a feature of the sample
            dist_function (function): A distance function which compute
                distance between two lines : dist_function(a,b) ->
                distance(a,b) where a and b are two lines of samples

        Returns:
            np.darray: The distance matrix
        """
        nb_samples = len(samples)
        distances = np.zeros((nb_samples, nb_samples))
        for i in range(nb_samples):
            for j in range(i):
                distances[i][j] = dist_function(samples[i], samples[j])
                distances[j][i] = distances[i][j]
        return distances

    def _radius_compatible_sets(self, distances, threshold):
        """Generate compatible sets from a distance matrix.

        Compatible sets are defined by the availability of a representative
        element. Each compatible set is associated with a set of centers
        stored in the same index as the compatible set in the attribute _centers.

        Args:
            distances (numpy.ndarray): A distance matrix
            threshold (float): The maximum distance at which elements are
                considered compatible

        Returns:
            list: A list of sets where there exist at least one elemen
                compatible with every other element
        """
        self._centers = []
        if len(distances) == 0:
            return []

        ind_line, ind_col = np.where(distances <= threshold)
        compatibles = []
        size = len(distances)
        for elem in range(size):
            compatible = frozenset(ind_col[ind_line == elem])
            is_in = False
            for i, comp_set in enumerate(compatibles):
                if comp_set == compatible:
                    self._centers[i].append(elem)
                    is_in = True
                    break
            if not is_in:
                compatibles.append(compatible)
                self._centers.append([elem])
        return compatibles

    def _set_membership(self, sets):
        """Compute set membership.

        Transform sets of elements into sets of memberships.
        Ouput set i contains element j if and only if input
        set j contains element i.

        Args:
            sets (list): A list of sets of compatible elements.
                Elements of the sets should be positive integers.
                They should cover all numbers between 0 and the
                maximum element

        Returns:
            list: The list of membership sets

        Examples:
            >>> _set_membership([{0, 1}, {0, 2}])
            [{0, 1}, {0}, {1}]
        """
        population = set().union(*sets)

        cover = [set() for _ in population]
        for index, set_ in enumerate(sets):
            for i in set_:
                cover[i].add(index)
        return cover

    def _greedy_min_set_cover(self, sets, candidate_centers=None):
        """Compute a minimal set cover.

        This implementation uses a greedy algorithm. The result is therefore
        not guaranteed to be minimal. Centers of compatible sets are kept.

        Args:
            sets (list): A list of compatible sets.
            candidate_centers (list): The list of centers associated with
                compatible sets. Default is None so attribute _centers will
                be used, no need to change.

        Returns:
            list: A list of one list of sets representing the minimal set cover.
        """
        greedy_cover = []
        unassigned = set().union(*sets)
        if candidate_centers is None:
            candidate_centers = self._centers
        centers = []
        candidates = list(sets)
        while len(unassigned) > 0:
            # find the largest covering set
            can_len = [len(set_ & unassigned) for set_ in candidates]
            ind_candidate = np.argmax(can_len)
            candidate = candidates.pop(ind_candidate)
            greedy_cover.append(frozenset(candidate))
            # or greedy_cover.append(frozenset(candidate & unassigned))
            unassigned = unassigned - candidate
            if self._measure == "radius":
                centers.append(candidate_centers.pop(ind_candidate))

        # Delete useless sets, starting by the smallest sets
        nb_sets = len(greedy_cover)
        if nb_sets > 1:
            for i in range(nb_sets - 1, -1, -1):
                other_elem = frozenset.union(
                    *(set(greedy_cover) - set({greedy_cover[i]}))
                )
                all_in = True
                for element in greedy_cover[i]:
                    if element not in other_elem:
                        all_in = False
                        break
                if all_in:
                    greedy_cover.pop(i)
                    if self._measure == "radius":
                        centers.pop(i)

        if self._measure == "radius":
            self._centers = [centers]
        return [greedy_cover]

    def _lp_min_set_cover(self, sets, cutoff=None, candidate_centers=None):
        """Compute a minimal set cover.

        This implementation uses linear programming. The result is guaranteed
        to be minimal if enough time is set. A cutoff can be set to improve
        computation time, but the result is therefore no guaranteed to be minimal.
        If cutoff is too low, solution may not be found at all.

        Args:
            sets (set): An iterable of sets.
            cutoff (float): Maximum runtime in seconds.
            candidate_centers (list): The list of centers associated with
                compatible sets. Default is None so attribute _centers will
                be used, no need to change.

        Returns:
            list: A list of one set of sets representing the minimal set cover.
        """
        if candidate_centers is None:
            candidate_centers = self._centers

        # Compute the population to cover
        population = set(frozenset.union(*sets))

        # determine all sets needed to cover population
        if self._solver is None:
            mod = Model("min set cover")
        else:
            mod = Model("min set cover", solver_name=self._solver)
        mod.verbose = 0

        # binary variables: s[i] is part of the solution if x_in[i]
        x_in = [mod.add_var(var_type=BINARY) for _ in range(len(sets))]
        # objective function: minimize the number of sets in the solution
        mod.objective = minimize(xsum(x_in[i] for i in range(len(sets))))
        # constraints: each element must be present in the solution
        for elem in population:
            mod += (
                xsum(x_in[i] for i, _set in enumerate(sets) if elem in _set) >= 1,
                str(elem) + " is covered",
            )

        # a cutoff can be set
        if cutoff is None:
            status = mod.optimize()
        else:
            status = mod.optimize(max_seconds=cutoff)
        if self._verbose == 1:
            if status == OptimizationStatus.OPTIMAL:
                print(
                    "Linear programming found an optimal min set cover with {} set(s).".format(
                        mod.objective_value
                    )
                )
            elif status == OptimizationStatus.FEASIBLE:
                print(
                    "Linear programming found a min set cover with {} set(s) but best possible with {}.".format(
                        mod.objective_value, mod.objective_bound
                    )
                )
            elif status == OptimizationStatus.NO_SOLUTION_FOUND:
                print(
                    "no feasible solution found, lower bound is: {}".format(
                        mod.objective_bound
                    )
                )

        covering = []
        centers = []
        for index, elem in enumerate(sets):
            if x_in[index].x:
                covering.append(elem)
                if self._measure == "radius":
                    centers.append(candidate_centers[index])

        if self._measure == "radius":
            self._centers = [centers]

        return [covering]

    def _select_candidate_minimize_distances(
        self, covers, distances, candidate_centers=None
    ):
        """Compute the final min set cover, i.e. the clustering, by choosing
        the best candidate : the one which minimize the intra-cluster distance
        and by affecting elements to only one cluster according to this objective.

        Args:
            covers (list): A list of different covers possible.
            distances (numpy.ndarray): A distance matrix.
            candidate_centers (list): The list of centers associated with
                compatible sets. Default is None so attribute _centers will
                be used, no need to change.

        Returns:
            set: The final clustering
        """
        # If covers is empty, return an empty set.
        if not covers:
            return set()

        if candidate_centers is None:
            candidate_centers = self._centers

        # For each possible cover :
        #   Find elements which are in multiple sets
        #   Compute the medoids
        #   Affect elements to the set with the closest medoid
        solutions = []
        for index_sol, cover in enumerate(covers):
            no_clusters_medoid = set()
            candidate_cluster = [set(e) for e in cover]
            membership = self._set_membership(candidate_cluster)
            # Delete elements that are in multiple sets
            for index, elem in enumerate(membership):
                if len(elem) > 1:
                    for no_sets in elem:
                        no_clusters_medoid.add(no_sets)
                        candidate_cluster[no_sets] = candidate_cluster[no_sets] - {
                            index
                        }

            # If measure is "radius", keep at least one center in each compatible set
            if self._measure == "radius":
                for index_cluster, centers in enumerate(candidate_centers[index_sol]):
                    # First, check if the compatible set contains at least one center
                    is_in = False
                    for center in centers:
                        if center in candidate_cluster[index_cluster]:
                            is_in = True
                            break
                    # If it does not contain any center, then add the one which minimize distances
                    if not is_in:
                        subdistance = distances[centers][:, centers]
                        ind_medoid = np.argmin(subdistance.sum(axis=1))
                        medoid = centers[ind_medoid]
                        candidate_cluster[index_cluster].add(medoid)
                        # Remove medoid from element that are in multiple sets
                        membership[medoid] = {index_cluster}

            # Compute medoid
            medoids = []
            for ind, subset in enumerate(candidate_cluster):
                # Only for sets that contained element in multiple sets
                if ind in no_clusters_medoid:
                    sublist = list(subset)
                    # Because we can look at the whole matrix
                    subdistance = distances[sublist][:, sublist]
                    ind_medoid = np.argmin(subdistance.sum(axis=1))
                    medoids.append(sublist[ind_medoid])
                else:
                    medoids.append(-1)

            # Affect elements
            for index, elem in enumerate(membership):
                if len(elem) > 1:
                    min_dist = -1
                    for no_sets in elem:
                        dist = distances[index, medoids[no_sets]]
                        if dist < min_dist or min_dist < 0:
                            min_dist = dist
                            no_cluster = no_sets
                    candidate_cluster[no_cluster] = candidate_cluster[no_cluster].union(
                        {index}
                    )

            solutions.append(set(frozenset(e) for e in candidate_cluster))

        # If there is only one possible solution, return it.
        if len(solutions) == 1:
            return solutions[0]

        # Else, find the best solution
        errors = []
        for sol in solutions:
            # Compute clustering error
            error = 0
            for subset in sol:
                # Compute medoid
                sublist = list(subset)
                subdistance = distances[sublist][:, sublist]
                ind_medoid = np.argmin(subdistance.sum(axis=1))
                medoid = sublist[ind_medoid]

                # Compute error
                for element in subset:
                    error += distances[element, medoid] ** 2
            errors.append(error)

        # Problem : Do not choose among equivalent solutions
        return solutions[np.argmin(errors)]

    def _compute_labels(self, clustering):
        """Create an array of cluster labels.

        Element i is in cluster j means that i-th element of the array is j.

        Args:
            clustering (set): A set of all clusters.
                Each element should be in only one cluster. Elements
                should be positive integers, and they should cover all
                numbers between 0 and the maximum element.

        Returns:
            list: A list of all cluster labels.
        """
        population = set().union(*clustering)
        cover = [None] * len(population)
        for index, set_ in enumerate(clustering):
            for i in set_:
                cover[i] = index
        return cover

    def clustering(self, input_mat, local_threshold=None, distance_function=None):
        """Compute a clustering with data represented as files and
        respecting the maximum intra-cluster distance set by the threshold.

        Args:
            input_mat (np.darray or list): Either a distance matrix, either
                an array of samples where each ligne is a sample and each
                column is a feature of the sample. If it is a distance matrix,
                distance_function should be None, else, distance function
                should be set. If it is a distance matrix, only the whole
                matrix will be read, so the distance should be symmetric
                and the distance between one point and itself should be zero.
            local_threshold (float): Distance threshold, i.e. the maximum
                diameter of clusters for the max-diameter algorithm or the
                maximum distance to the cluster representative for the max-radius
                algorithm. If threshold has not been set in constructor,
                this parameter is needed. Else, it is optional, it won't
                overwrite the threshold of the constructor but local_threshold
                will be use instead if set. Defaults to None.
            distance_function (function): A distance function to compute a
                distance matrix from input_mat. If None, it means that input_mat
                is already a distance matrix.

        Returns:
            list: A list of cluster labels.

        Raises:
            TypeError: When threshold is not a number (integer or float).
            ValueError: When attribute measure is not "diameter" or "radius" or
                cover is not "greedy", "exact" or "lp".
        """
        # Choose and check the threshold
        if local_threshold is not None:
            if type(local_threshold) in {int, float}:
                threshold = local_threshold
            else:
                raise TypeError("local_threshold should be a number (integer or float)")
        else:
            if self._threshold is not None:
                threshold = self._threshold
            else:
                raise ValueError(
                    "local_threshold is needed because threshold attribute is None"
                )

        if self._verbose == 1:
            print(f"Run clustering with a {self._measure} constraint of {threshold}.")
            print(f"Minimum set cover is computed using {self._cover} algorithm.")

        # Read or compute distance matrix
        if distance_function is None:
            # Make sure that distance matrix it is a numpy array
            distances = np.array(input_mat)
        else:
            distances = self._compute_distances(input_mat, distance_function)

        # Compute compatible sets
        if self._measure == "radius":
            compatible_sets = self._radius_compatible_sets(distances, threshold)
        else:
            raise ValueError("measure should be either radius or diameter")

        if self._verbose == 2:
            print("nb_compatible_sets:" + str(len(compatible_sets)))

        # Compute min set cover
        if self._cover == "greedy":
            covering_sets = self._greedy_min_set_cover(compatible_sets)
        elif self._cover == "lp":
            covering_sets = self._lp_min_set_cover(compatible_sets, self._cutoff)
        else:
            raise ValueError("cover should be greedy, exact or lp")

        if self._verbose == 2:
            print("nb_solutions:" + str(len(covering_sets)))

        # Compute final clustering
        clustering = self._select_candidate_minimize_distances(
            covering_sets, distances
        )

        if self._verbose == 2:
            print("nb_clusters:" + str(len(clustering)))

        if self._verbose == 1:
            print(f"Final clustering contains : {len(clustering)} cluster(s).")

        # Compute cluster labels
        labels = self._compute_labels(clustering)

        return labels
