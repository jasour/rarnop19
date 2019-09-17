% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 3: Sum Of Squares Based SDP For Nonlinear Optimization
%% SOS based SDP for Unconstrained optimization : min_x p(x)

clc;clear all
% variable x and lower bound gamma
sdpvar x  gamma
% objective function p(x)
p = x^4+2*x^3-12*x^2-2*x+6;
% SOS condition
F = sos(p-gamma);
% SDP Solver
ops = sdpsettings('solver','mosek');
% Solve SOS based SDP
[sol,v,Q]=solvesos(F,-gamma,ops,gamma);

%% Obtained Result

% Obtained lower bound
value(gamma)
% Obtaied PSD matrix in SOS representation p(x)-gamma=B'(x)QB(x)
value(Q{1})
% obtained basis B(x)
sdisplay(v{1})
% obtains SOS representation: p(x)=h(x)'h(x)
h=sosd(F); sdisplay(h);sdisplay(h'*h)
% test p(x)-h(x)'h(x)=gamma
clean(p-h'*h,1e-6)