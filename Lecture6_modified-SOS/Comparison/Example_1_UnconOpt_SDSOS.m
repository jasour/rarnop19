% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 6: Modified Sum-of-Squares Relaxations for Large Scale Optimizations
%%  SDSOS.DSOS based SDP for Unconstrained optimization : min_x p(x)

clc;clear all

% n: dimension of polynomial,  d: order of polynomial
n=30;d=2;

% polynomial variables
x = msspoly('x',n);
prog = spotsosprog;  prog = prog.withIndeterminate(x);

% objective function
p=5; for i=2:n;  p=p+100*(x(i)-x(i-1).^2).^2+(1-x(i)).^d; end

%lower bound gamma
[prog,gamma] = prog.newFree(1);

% SDSOS/DSOS condition
prog = prog . withSDSOS (p-gamma) ;

tic
sol = prog . minimize ( -gamma,@spot_mosek);
t=toc
%% Obtained Result
% Obtained lower bound
Optimal_Objective=double(sol.eval(gamma))
Time=toc

