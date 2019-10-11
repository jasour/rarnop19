% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 3: Sum Of Squares Based SDP For Nonlinear Optimization
%% SOS Decomposition
% Dx/dt=f*x, x=[x1;x2;x3]
% f=[(-x1^3-x1*x3^2);(-x2-x1^2*x2);(-x3-3*x3/(x3^2+1)+3*x1^2*x3)]
% V(x)= 5*x1^2+4*x2^2+x3^2;
% Show that nonlinear systemis is stable.

clc;clear
% variables
 sdpvar x1 x2 x3

% V(x)
V = 5*x1^2+4*x2^2+x3^2;

%f=[(-x1^3-x1*x3^2);(-x2-x1^2*x2);(-x3-3*x3/(x3^2+1)+3*x1^2*x3)];
% f2=f1*(x3^2+1)
f2=[(-x1^3-x1*x3^2)*(x3^2+1);(-x2-x1^2*x2)*(x3^2+1);-x3*(x3^2+1)-3*x3+(3*x1^2*x3)*(x3^2+1)];


x = [x1;x2;x3];
% -dv/dt*((x3^2+1))
Dv =-jacobian(V,x)*f2;

%% SOS SDP
F = sos(Dv);
[sol,v,Q]=solvesos(F);

%% p(x)=B(x)'Q B(x)
% obtained PSD matrix Q: p(x)=B(x)'Q B(x)
Q=value(Q{1});
% obtained basis B(x)
sdisplay(v{1});
% obtains sum of squares polynomials: p(x)=h(x)'h(x)
h=sosd(F); sdisplay(h)
% test p(x)-h(x)'h(x)=0
clean(Dv-h'*h,1e-6)

% Hence, V(x) is a Lyapunov function (V(0)=0, V(x)>=0, -dV(x)/dt) <=0 ),
% and system is stable.

