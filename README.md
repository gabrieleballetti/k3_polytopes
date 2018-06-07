# K3 polytopes

[K3 Polytopes and their Quartic Surfaces](https://arxiv.org/pdf/1806.02236.pdf) - Gabriele Balletti, Marta Panizzut, Bernd Sturmfels

Source code for the classification of K3 polytopes. The code is written for [polymake](https://polymake.org/), and uses its implementation of [TOPCOM](http://www.rambau.wm.uni-bayreuth.de/TOPCOM/).

## Content

* find_subpolytopes.pl - Find all the subpolytopes of the fourth dilation of the standard simplex up to permutation of the coordinates
* find_triangulations.pl - Extract all the central regular unimodular triangulations of the input polytope
* generate_k3.pl - Given a reflexive polytope and a central regular unimodular triangulation, computes the K3 polytope defined by the qurtic surface dual to the triangulation.
* count_combinatorial_types.pl - Extracts the different vertex-facets incidence graphs of K3 polytopes
