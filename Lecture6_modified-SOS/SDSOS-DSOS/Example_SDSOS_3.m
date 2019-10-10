% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 3: Sum Of Squares Based SDP For Nonlinear Optimization
%% Scaled Diagonally SOS (SDSOS) and Diagonally  SOS (DSOS) based SDP for Constrained optimization : min_x p(x), s.t. g_i(x)>=0 i=1,...,m

clc;clear all

% Relaxation Order
d=3;

% variable x 
x = msspoly('x',2);

% DSOS/SDSOS Programing
prog = spotsosprog;  prog = prog.withIndeterminate(x);

% objective function
p=x(1)^6/3 - (21*x(1)^4)/10 + 4*x(1)^2 + x(1)*x(2) + 4*x(2)^4 - 4*x(2)^2 + 3/2 ;
% Constraints g_i(x)>=0, i=1,2,3
g=-x(1)^4/16 + x(1)^3/4 - x(1)^2/4 - (9*x(2)^2)/100 + 29/400;

%lower bound gamma
[prog,gamma] = prog.newFree(1);

% vector of monomials up to order d
mos=monomials(x,0:2*d); 
%sigma_1
[prog,coeffs1] = prog.newFree(length(mos)); s1 = coeffs1'*mos;

%dsos/sdsos/sos conditions: sos(p-gamma-[s1 s2 s3]*g), sos(s0), sos(s1), sos(s2), sos(s3);
prog = prog . withDSOS (p-gamma-[s1]*g) ;
prog = prog . withDSOS (s1) ;


% SOlve SDP  SDP solvers: @spot_sedumi, @spot_mosek
sol = prog . minimize ( -gamma,@spot_sedumi) ;

% obtaind gamma
double(sol.eval(gamma))

