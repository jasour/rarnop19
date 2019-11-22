% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 6: Modified Sum-of-Squares Relaxations for Large Scale Optimizations
%% Sparse SOS based SDP for Sparse Constrained optimization 

clc;clear all;


% n: dimension of polynomial,  d: order of polynomial
n=10;d=6;

%%
% variable x and lower bound gamma
x = sdpvar(n,1); sdpvar gamma

% objective function p(x)=5+ Sum_{i=1}^n {  100(x(i+1)-x(i))^2 + (x(i)-1)^2 }
p=5; for i=1:n-1;  p=p+100*(x(i+1)-x(i).^2).^2+(1-x(i)).^d; end

% constraints -u<=x(i)<=u i=1,...,N
g=[2-x.^2];

% relaxation order
dr = 1; 

% Sparse SOS polynomials sigma_i
C=[];%coefficient vector of s polynomials
for i=1:n;  
    % s(i)  is Sigma_i(x) polynomial of order d with Coeffs c
    [s(i),c(:,i)] = polynomial(x(i),2*dr);  C=[C;c(:,i)];
end

% SOS Conditions:
F = [sos(p - gamma -[s]*g), sos(s),[C;gamma]];

% SDP Solver
ops = sdpsettings('solver','mosek');

% Chordal Sparsity 
ops.sos.csp = 1;   

% Solve SDP
[sol,v,Q]   = solvesos(F,-gamma,ops);

%% Obtained Result
% Obtained lower bound
Optimal_Objective=value(gamma)
Time=sol.solvertime   






