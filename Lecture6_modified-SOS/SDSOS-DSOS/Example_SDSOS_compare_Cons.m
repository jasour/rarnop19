% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 3: Sum Of Squares Based SDP For Nonlinear Optimization
%% SOS based SDP for Constrained optimization 

% number of the variables
N=10;

% Relaxation Order
d=1;

x = msspoly('x',N);
prog = spotsosprog;  prog = prog.withIndeterminate(x);


p=5; g=0*x;
for i=1:N
   p=p+(x(i)-1).^2; % objective function% objective function
   g(i)=1-x(i)^2; % constraint
end

%gamma
[prog,gamma] = prog.newFree(1);

mos=monomials(x,0:2*d);
coeffs=gamma*zeros(length(mos),N);s=gamma*zeros(N,1);
for i=1:N  
[prog,coeffs(:,i)] = prog.newFree(length(mos));
s(i,1) = coeffs(:,i)'*mos;
end

%F = [sos(p-gamma-s0-sum{si(x)*gi(x)}], sos(si(x)) i=1,...,m];
prog = prog . withSDSOS (p-gamma-[s']*g) ;
prog = prog . withSDSOS (s) ;

tic
sol = prog . minimize ( -gamma,@spot_mosek);
t=toc
double(sol.eval(gamma))

