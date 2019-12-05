# Multi-scan trajectory PMBM filter
MATLAB implementations of the Multi-Scan Trajectory Poisson Multi-Bernoulli Mixture Filter, Multi-Bernoulli Mixture Filter and Multi-Bernoulli Mixture-01 Filter

This repository contains the Matlab implementations of the Multi-Scan Trajectory Poisson Multi-Bernoulli Mixture Filter via dual decomposition proposed in 

Y. Xia, K. Granström, L. Svensson and Á. F. García-Fernández, "An Implementation of the Poisson Multi-Bernoulli Mixture Trajectory Filter via Dual Decomposition," 2018 21st International Conference on Information Fusion (FUSION), Cambridge, 2018

Full text is available at https://arxiv.org/pdf/1811.12281.pdf

and in

Y. Xia, K. Granström, L. Svensson, Á. F. García-Fernández and Jason L. Williams, "Multi-Scan Implementation of the Trajectory Poisson Multi-Bernoulli Mixture Filter," Accepted for publication in Journal of Advances in Information Fusion (special issue on Multiple Hypothesis Tracker), 2020

Full text is available at https://arxiv.org/abs/1912.01748

More details on the point target PMBM tracker can be found in the paper 

Granström, Karl, Lennart Svensson, Yuxuan Xia, Jason Williams, and Ángel F. García-Femández. "Poisson multi-Bernoulli mixture trackers: continuity through random finite sets of trajectories." In 2018 21st International Conference on Information Fusion (FUSION). IEEE, 2018.

Full text is available at https://arxiv.org/pdf/1812.05131.pdf

More details on Dual Decomposition can be found in the paper

N. Komodakis, N. Paragios, and G. Tziritas, "MRF energy minimization and beyond via dual decomposition," IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 33, no. 3, pp. 531–552, 2011.

The filters are evaluated using 

1. the generalised optimal subpattern-assignment (GOSPA)

Abu Sajana Rahmathullah, Ángel F. García-Fernández, Lennart Svensson, Generalized optimal sub-pattern assignment metric, in 20th International
Conference on Information Fusion, 2017.

Video on GOSPA: https://www.youtube.com/watch?v=M79GTTytvCM

2. the trajectory metric

Abu Sajana Rahmathullah, Ángel F. García-Fernández, Lennart Svensson, A metric on the space of finite sets of trajectories for evaluation of multi-target tracking algorithms

