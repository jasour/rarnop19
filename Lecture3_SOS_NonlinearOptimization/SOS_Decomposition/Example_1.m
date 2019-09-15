% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 3: Sum Of Squares Based SDP For Nonlinear Optimization
%% SOS Decomposition

clc;clear
% variables
x = sdpvar(1,2);
% polynomial p(x)
p = x(1)^2-x(1)*x(2)^2+x(2)^4+1;

%% SOS SDP
F = sos(p);
[sol,v,Q]=solvesos(F);

%% p(x)=B(x)'Q B(x)
% obtained PSD matrix Q: p(x)=B(x)'Q B(x)
Q=value(Q{1});
% obtained basis B(x)
sdisplay(v{1});
% obtains sum of squares polynomials: p(x)=h(x)'h(x)
h=sosd(F); sdisplay(h);
% test p(x)-h(x)'h(x)=0
clean(p-h'*h,1e-6)

%% Q=LL'
% Eigenvalue Decomposition  Q= V D V'
[V,D]=eig(Q);
% L = V*sqrt(D);
L=(V*sqrt(D));
% test Q-LL'=0
clean(Q-L*L',1e-6)
