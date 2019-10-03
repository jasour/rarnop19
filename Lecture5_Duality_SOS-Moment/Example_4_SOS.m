% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 3: Sum Of Squares Based SDP For Nonlinear Optimization
%% SOS based SDP for Constrained optimization : min_x p(x), s.t. g_i(x)>=0 i=1,...,m
clc;clear all

% variable x1, x2 and lower bound gamma
sdpvar x1 x2 gamma
% objective function p(x)
p=-x1^2*x2^2+x1^4*x2^2+x1^2*x2^4;
% added constraint
g=2-x1^2-x2^2;

% SOS d order polynomials sigma_i i=1,2,3 
% si: sigma_i  ci: coefficient vector, i=1,2,3
d=3;
[s1,c1] = polynomial([x1 x2],2*d-2);


% SDP Solver
ops = sdpsettings('solver','mosek');
% SOS Conditions
F = [sos(p-gamma-[s1]*g), sos(s1)];
% Solve SOS based SDP
[sol,v,Q]=solvesos(F,-gamma,ops,[c1;gamma]);

%% %% Obtained Results

% Obtained lower bound
value(gamma)

% Obtained coefficient vector ci: , i=1,2,3 for SOS si: sigma_i  
[value(c1)'];


%S1: SOS sigma_1  
S1=sosd(F(2))'*sosd(F(2));
sdisplay(v{2}'*Q{2}*v{2})
sdisplay(sosd(F(2))'*sosd(F(2)))

%SOS decomposition of p(x): h'h or B'(x)QB(x)
sdisplay(v{1}'*Q{1}*v{1})
h=sosd(F(1));
sdisplay(sosd(F(1))'*sosd(F(1)))

% test: p(x)-S1g_1(x)-S2g_2(x)-S3g_3(x) - h'h=gamma
clean(p-[S1]*g - h'*h,1e-5)
