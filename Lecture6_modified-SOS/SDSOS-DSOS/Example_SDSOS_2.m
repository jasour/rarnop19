% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 3: Sum Of Squares Based SDP For Nonlinear Optimization
%% Scaled Diagonally SOS (SDSOS) and Diagonally  SOS (DSOS) based SDP for Constrained optimization : min_x p(x), s.t. g_i(x)>=0 i=1,...,m

clc;clear all

% Relaxation Order
d=2;

% variable x 
x = msspoly('x',2);

% DSOS/SDSOS Programing
prog = spotsosprog;  prog = prog.withIndeterminate(x);

% objective function
p = (1+x(1)*x(2))^2-x(1)*x(2)+(1-x(2))^2 ;
%-x(1);
% Constraints g(x)>=0
g=[3-2*x(2)-x(1)^2-x(2)^2;-x(1)-x(2)-x(1)*x(2);1+x(1)*x(2)];

%lower bound gamma
[prog,gamma] = prog.newFree(1);

% vector of monomials up to order d
mos=monomials(x,0:2*d); 
%sigma_1
[prog,coeffs1] = prog.newFree(length(mos)); s1 = coeffs1'*mos;
%sigma_2
[prog,coeffs2] = prog.newFree(length(mos)); s2 = coeffs2'*mos;
%sigma_3
[prog,coeffs3] = prog.newFree(length(mos)); s3 = coeffs3'*mos;

%dsos/sdsos/sos conditions: sos(p-gamma-[s1 s2 s3]*g), sos(s0), sos(s1), sos(s2), sos(s3);
prog = prog . withDSOS (p-gamma-[s1 s2 s3]*g) ;
prog = prog . withDSOS (s1) ;
prog = prog . withDSOS (s2) ;
prog = prog . withDSOS (s3) ;

% SOlve SDP  SDP solvers: @spot_sedumi, @spot_mosek
sol = prog . minimize ( -gamma,@spot_mosek) ;

% obtaind gamma
double(sol.eval(gamma))

