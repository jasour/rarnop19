% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 3: Sum Of Squares Based SDP For Nonlinear Optimization
%% SOS based SDP for Constrained optimization : min_x p(x), s.t. g_i(x)>=0 i=1,...,m
clc;clear all

% variable x1, x2 and lower bound gamma
sdpvar x1 x2 gamma
% objective function p(x)
p = -x1;
% Constraints g_i(x)>=0, i=1,2,3
g=[3-2*x2-x1^2-x2^2;-x1-x2-x1*x2;1+x1*x2];

% SOS d order polynomials sigma_i i=1,2,3 
% si: sigma_i  ci: coefficient vector, i=1,2,3
d=2;
[s1,c1] = polynomial([x1 x2],2*d-2);
[s2,c2] = polynomial([x1 x2],2*d-2);
[s3,c3] = polynomial([x1 x2],2*d-1);

% SDP Solver
ops = sdpsettings('solver','mosek');
% SOS Conditions
F = [sos(p-gamma-[s1 s2 s3]*g), sos(s1), sos(s2), sos(s3)];
% Solve SOS based SDP
[sol,v,Q]=solvesos(F,-gamma,ops,[c1;c2;c3;gamma]);

%% %% Obtained Results

% Obtained lower bound
P_opt=value(gamma)
Qopt=Q{1};
% Obtained coefficient vector ci: , i=1,2,3 for SOS si: sigma_i  
%[value(c1)';value(c2)';value(c3)'];


%S1: SOS sigma_1  
S1=sosd(F(2))'*sosd(F(2));
sdisplay(v{2}'*Q{2}*v{2})
sdisplay(sosd(F(2))'*sosd(F(2)))

%S2: SOS sigma_2  
S2=sosd(F(3))'*sosd(F(3));
sdisplay(v{3}'*Q{3}*v{3})
sdisplay(sosd(F(3))'*sosd(F(3)))

%S3: SOS sigma_3  
S3=sosd(F(4))'*sosd(F(4));
sdisplay(v{4}'*Q{4}*v{4})
sdisplay(sosd(F(4))'*sosd(F(4)))

%SOS decomposition of p(x): h'h or B'(x)QB(x)
sdisplay(v{1}'*Q{1}*v{1})
h=sosd(F(1));
sdisplay(sosd(F(1))'*sosd(F(1)))

% test: p(x)-S1g_1(x)-S2g_2(x)-S3g_3(x) - h'h=gamma
clean(p-[S1 S2 S3]*g - h'*h,1e-5)
