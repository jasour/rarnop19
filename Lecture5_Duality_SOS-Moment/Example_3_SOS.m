% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 3: Sum Of Squares Based SDP For Nonlinear Optimization
%% SOS based SDP for Unconstrained optimization : min_x p(x)

clc;clear all
% variable x and lower bound gamma
sdpvar x1 x2  gamma
% objective function p(x)
p = -x1^2*x2^2+x1^4*x2^2+x1^2*x2^4;
% SOS condition
F = sos(p-gamma);
% SDP Solver
ops = sdpsettings('solver','sedumi');
% Solve SOS based SDP
[sol,v,Q]=solvesos(F,-gamma,ops,gamma);

%% Obtained Result

% Obtained lower bound
value(gamma)
% Obtaied PSD matrix in SOS representation p(x)-gamma=B'(x)QB(x)
value(Q{1});
% obtained basis B(x)
sdisplay(v{1});
% obtains SOS representation: p(x)=h(x)'h(x)
h=sosd(F); sdisplay(h);sdisplay(h'*h);
% test p(x)-h(x)'h(x)=gamma
clean(p-h'*h,1e-6);

%% (p(x)-gamma_s) does not admit SOS representation
F = sos(p-1/27);solvesos(F);
h = sosd(F); sdisplay(h)