% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 3: Sum Of Squares Based SDP For Nonlinear Optimization
%% Sparse SOS based SDP for Sparse Constrained optimization 

clear all;clc;close all

N = 200; % number of variables
% variable x and lower bound gamma
x = sdpvar(N,1); sdpvar gamma

% objective function p(x)=5+ Sum_{i=1}^n {x(i)-1)^2}
p=5; for i=1:N;  p=p+(x(i)-1)^2; end


% SOS Conditions:
F = [sos(p - gamma), [gamma]];

% SDP Solver
ops = sdpsettings('solver','mosek');

% Chordal Sparsity 
ops.sos.csp = 1;   

% Solve SDP
[sol,v,Q]   = solvesos(F,-gamma,ops);

% Obtained Results
t_csp       = sol.solvertime   
lower   = value(gamma)





