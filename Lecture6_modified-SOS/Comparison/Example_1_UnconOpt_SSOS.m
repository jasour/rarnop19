% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 6: Modified Sum-of-Squares Relaxations for Large Scale Optimizations
%% Sparse SOS based SDP for Unconstrained optimization : min_x p(x)

clc;clear all
% n: dimension of polynomial,  d: order of polynomial
n=10;d=6;

%%
% polynomial variables
x=sdpvar(1,n);
% lower bound
sdpvar  gamma
% objective function p(x)
p=5; for i=2:n;  p=p+100*(x(i)-x(i-1).^2).^2+(1-x(i)).^d; end

% SOS condition
F = sos(p-gamma);
% SDP Solver
ops = sdpsettings('solver','mosek');

% cordal sparsity
ops.sos.csp = 1;  

% Solve SOS based SDP
tic
[sol,v,Q]=solvesos(F,-gamma,ops,gamma);
toc
%% Obtained Result
% Obtained lower bound
Optimal_Objective=value(gamma)
Time=toc