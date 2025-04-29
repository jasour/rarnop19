# Julia Code for Chance Constrained Sets


**Chance Constrained Set:** Set of all design parameters (control inputs, plans,...) that satisfy probabilistic nonlinear safety constraints.

- By constructing chance constrained sets, we can turn a chance constrained optimization in to a deterministic optimization.
- Moment-sum-of-squates based Semidefinite Programming
- More Information: Lecture 7: Ashkan Jasour, "Risk Aware and Robust Nonlinear Planning", MIT 16.S498, 2019.
- https://rarnop.mit.edu/Lectures-Codes
- https://github.com/jasour/Chance-Optimization



- CCset_Julia_MomentOpt_1.jl: Example 1: Lecture 7, Page 133



- To run the codes you need the following packages in Julia: 1) MomentOpt julia package : Pkg.add("MomentOpt"), 2) DynamicPolynomials: Pkg.add("DynamicPolynomials"), 3) SDP solvers like "mosek": Pkg.add("MosekTools")

- To plot the solution of dual SDP run Main_Plot.m
