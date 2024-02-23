# Equiwide Clustering

The **Equiwide Clustering** algorithm computes a clustering based on a maximum intra-cluster distance. In this algorithm, two elements can only be assigned to the same cluster if they are separated by a distance of less than a given threshold. This threshold is called the maximum diameter. A variant of the algorithm considers a maximum radius instead of a maximum diameter: each cluster then has a medoid whose distance from each element of the same cluster is lower than the provided threshold.
The main objective is to minimize the number of clusters under this constraint.

As there can be many solutions providing the same number of clusters, the selected clustering is determined by either minimizing the intra-cluster distance or by maximizing the size of the largest clusters.

This algorithm forms clusters with bounded width, i.e. these clusters define equi-wide regions.

This algorithm is based on the following paper: [Clustering to the Fewest Clusters Under Intra-Cluster Dissimilarity Constraints](https://hal.science/hal-03356000).

## Contributors

- [Jennie Andersen](https://liris.cnrs.fr/page-membre/jennie-andersen)
- [Brice Chardin](https://www.lias-lab.fr/members/bricechardin/)
- [Mickael Baron](https://www.lias-lab.fr/members/mickaelbaron/)