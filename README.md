# K3 polytopes

Source code for the classification of K3 polytopes for the paper [K3 Polytopes and their Quartic Surfaces](https://arxiv.org/pdf/1806.02236.pdf). The code is written for [polymake](https://polymake.org/), and uses its implementation of [TOPCOM](http://www.rambau.wm.uni-bayreuth.de/TOPCOM/).

## Content

* find_subpolytopes.pl - find all the subpolytopes of the fourth dilation of the standard simplex up to permutation of the coordinates
* find_triangulations.pl - extract all the central regular unimodular triangulations of the input polytope
* generate_k3.pl - given a reflexive polytope and a central regular unimodular triangulation, computes the K3 polytope defined by the qurtic surface dual to the triangulation.
* count_combinatorial_types.pl - extracts the different vertex-facets incidence graphs of K3 polytopes
* data/ - folder where all the classifications and data are stored
