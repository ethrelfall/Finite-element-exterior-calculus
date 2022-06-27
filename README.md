# Finite-element-exterior-calculus
Repository for FEEC material e.g. Firedrake examples.

# Finite element exterior calculus in Firedrake

These are some simple Firedrake examples of the use of finite element exterior calculus as presented in the textbook Finite Element Exterior Calculus by D.N. Arnold.

All current examples are cases of the vector Laplacian on a planar two-dimensional domain.

In all cases, the .msh file has the same filename as the script.

**re-entrant_corner.py** - Laplacian problem with source on a L-shaped domain containing a 90 degree re-entrant corner.  Theory indicates in the continuous solution a r^{-1/3} singularity at distance r from the corner and the numerical challenge is to resolve this singularity without spurious oscillations.

**source_problem_on_annulus.py** - Laplacian problem with source on an annular domain.  The numerical challenge is to avoid pollution of the output by the 1-cohomology of the vector Laplacian on the domain (represented by the usual (-y/r^2, -x/r^2) vortex solution).

**eigenvalue_problem_single_vortex.py** - eigenvalue problem to find the 1-cohomology of the vector Laplacian on the annulus (same domain as source_problem_on_annulus.py script and the zero-eigenvalue mode is the vortex described for that example).

**eigenvalue_problem_double_vortex.py** - slightly more interesting eigenvalue problem to find the 1-cohomology of the vector Laplacian on a square with two holes removed.  The two zero-eigenvalue modes correspond roughly (but not precisely) to the cases of a co-rotating and a counter-rotating vortex pair (the solver does not find the pure co-/counter-rotating modes but rather linear combinations thereof). 

The help of Colin Cotter, Patrick Farrell, and Lawrence Mitchell is acknowledged.  Any mistakes are due to the author (Ed Threlfall). 
