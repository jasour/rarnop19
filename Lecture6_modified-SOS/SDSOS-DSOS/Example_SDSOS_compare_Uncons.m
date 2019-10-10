% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 3: Sum Of Squares Based SDP For Nonlinear Optimization
%% SOS based SDP for Unconstrained optimization 
clc;clear;close all

% number of the variables
N=200;

x = msspoly('x',N);
prog = spotsosprog;  prog = prog.withIndeterminate(x);


p=5; g=0*x;
for i=1:N
   p=p+(x(i)-1).^2; % objective function% objective function
end

%gamma
[prog,gamma] = prog.newFree(1);

% p-gamma in DSOS/SDSOS/SOS
prog = prog . withDSOS (p-gamma) ;


tic
% SOlve SDP
sol = prog . minimize ( -gamma,@spot_mosek);
t=toc

% obtaind gamma
double(sol.eval(gamma))

