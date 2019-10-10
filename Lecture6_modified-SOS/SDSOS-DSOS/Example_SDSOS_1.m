% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 3: Sum Of Squares Based SDP For Nonlinear Optimization
%% Scaled Diagonally SOS (SDSOS) and Diagonally  SOS (DSOS) based SDP for Unconstrained optimization : min_x p(x)

clc;clear all

% variable x 
x = msspoly('x',2);

% DSOS/SDSOS Programing
prog = spotsosprog;  
prog = prog.withIndeterminate(x);

% objective function
p = 3+2*x(1)+2*x(2)+3*x(1)^2+2*x(1)*x(2)+3*x(2)^2+x(1)^4+x(2)^4;

%lower bound gamma
[prog,gamma] = prog.newFree(1);

% p-gamma in DSOS/SDSOS/SOS
prog = prog . withDSOS (p-gamma) ;

% SOlve SDP
sol = prog . minimize ( -gamma,@spot_mosek) ;

% obtaind gamma
double(sol.eval(gamma))


