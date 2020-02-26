## A parallelized inverse fast multipole method (IFMM)

** Overview **

The code provides an (approximate) direct solver for the system of equations Ax=b, where the matrix entries a_{ij} are based on a kernel K(r), i.e. a_{ij} = K(||x_i - x_j||). The algorithm has a linear scaling O(N) and is based on two key ideas: (1) the original dense matrix A is transformed into an extended sparse matrix, and (2) fill-ins arising during the elimination phase are compressed as low-rank matrices if they correspond to well-separated interactions (hence ensuring that the sparsity pattern of the extended matrix is preserved). The method can be used as a highly accurate direct solver (machine accuracy can be achieved if needed), but it is more efficient to use it as a low accuracy preconditioner in an iterative solver such as GMRES.

---

** Description of the main code **

Please have a look at file main.cpp, which is an example of how to use the method. In this example, N points are randomly distributed in the unit cube. The most important aspects of this file are summarized below:

* Lines 14-16: input parameters (see below).

* Line 25: distribution of the points in a 3D cube.

* Line 28: defines the kernel K(r).

* Lines 30-34: iFMM parameters (see below).

* Lines 40-55: For validation purposes, we choose a vector 'x_exact' and calculate the corresponding right-hand side (b=Ax_exact). (We don't need to do this in an actual problem.)

* Lines 71-102: Solve Ax=b with the inverse multipole method (and compare this 'x' with 'x_exact'). As you can see on line 77, I create an object 'ifmm' of the class 'IFMM_Matrix'. What follows are initialization, setRHS, elimination, and substitution. Most of the time is spent in the elimination phase (certainly for large systems).

---

** Edit Makefile **

We assume to use g++, but the Intel's C++ compiler also can work.

* Specify the directory of the Eigen by EIGEN_INCLUDE.
* Eigen works better with MKL. To do so, add EIGEN_USE_MKL_ALL to FLAGS.
* Use the OMP library by -fopenmp.
* To parallelize the code, add "PARA_CHECK_RANK" and "IFMM_PARALLELIZE" to FLAGS.
* To parallelize the code, add "PARA_CHECK_RANK" and "IFMM_PARALLELIZE" to FLAGS.

---

** Parameters for IFMM **

Once the elimination has been performed, the factorization can easily (and fast) be re-used to solve to system for another right-hand side. This is illustrated on lines 105-123 (for 10 random RHS vectors): you need to reset the RHS-vector, re-use the factorization, and perform the substitution once more. Re-using the factorization is particularly useful if the iFMM is applied as a preconditioner in an iterative solver (which is not shown in this example).

As you will notice, some parameters have to be specified for the method:

* n: the number of Chebyshev points that are used in each direction for constructing low-rank approximations. The larger n, the more accurate the solver is, but at a higher computational cost. I would suggest using n=1, 2, or 3 for preconditioning purposes.
* l: the number of levels in the octree (must be larger than two). This will depend on the total number of points N.
* epsilon_rel_fillin: relative accuracy for compressing the fill-ins that arise during the elimination. You can play with this value - I suggest something between 10^{-1} and 10^{-3}.
* epsilon_rel_basis: relative accuracy for updating the basis of the FMM operators. 10^{-3} seems to be a good choice.
* LR_mode: choose LR_mode=3 for best performance (this uses a randomized singular value decomposition for obtaining low rank approximations).
    Delay_update: currently not supported yet, so put false (I still have some work to do...)

---

** References **

* Toru Takahashi, Chao Chen, Eric Darve, "Parallelization of the inverse fast multipole method with an application to boundary element method", Computer Physics Communications, 247, 2019. https://www.sciencedirect.com/science/article/pii/S0010465519303194

* Toru Takahashi, Pieter Coulier, Eric Darve, "Application of the inverse fast multipole method as a preconditioner in a 3D Helmholtz boundary element method", Journal of Computational Physics, 34, pp.406-428. https://www.sciencedirect.com/science/article/pii/S0021999117302875

* Pieter Coulier, Hadi Pouransari, Eric Darve, "The Inverse Fast Multipole Method: Using a Fast Approximate Direct Solver as a Preconditioner for Dense Linear Systems", SIAM J. Sci. Comput., 39(3), A761–A796.

* S. Ambikasaran, E. Darve, "The inverse fast multipole method", arXiv:1407.1572
, 2014. http://arxiv.org/abs/1407.1572

---

** Authors **

Toru Takahashi
Chao Chen
Pieter Coulier
Eric Darve
