% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 3: Sum Of Squares Based SDP For Nonlinear Optimization
%% Sparse SOS based SDP for Sparse Constrained optimization 

clear all;clc;close all

N = 20; % number of variables
% variable x and lower bound gamma
x = sdpvar(N,1); sdpvar gamma

% objective function p(x)=5+ Sum_{i=1}^n {x(i)-1)^2}
p=5; for i=1:N-1;  p=p+(x(i)-1).^2; end

% constraints -u<=x(i)<=u i=1,...,N
g=[2-x.^2];


% relaxation order
d = 1; 
%% Sparse SOS polynomials sigma_i
C=[];%coefficient vector of s polynomials
for i=1:N;  
    % s(i)  is Sigma_i(x) polynomial of order 2d with Coeffs c in terms of
    % variables of g_i(x)
    [s(i),c(:,i)] = polynomial(x(i),2*d);  C=[C;c(:,i)];
end


% SOS Conditions:
%(p-gamma- Sum_i Sigma_i(x)g_i(x)) \in SOS   and
% Sigma_i(x) \i n SOS
F = [sos(p - gamma -[s]*g), sos(s),[C;gamma]];

% SDP Solver
ops = sdpsettings('solver','mosek');

% Chordal Sparsity 
ops.sos.csp = 1;   

% Solve SDP
[sol,v,Q]   = solvesos(F,-gamma,ops);

% Obtained Results
t_csp       = sol.solvertime   
lower_csp   = value(gamma)





