% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 2: Nonlinear Optimization

%% 

clc;clear;close all;


opts = optimset('fmincon');opts.Display = 'iter';opts.Algorithm = 'active-set';opts.MaxFunEvals = 2000;



fun = @(x)(4-2.1*x(1).^2+(x(1).^4)/3).*x(1).^2 + x(1).*x(2) + (-4+4*x(2).^2).*x(2).^2

lb=[-2 -2]; ub=[2 2];

x0=[-1 -1];% x1=-0.0899  x2=0.7126
%x0=[1 1];% x1=0.0899 ,  x2=-0.7126
%x0=[1.5 -1]; % x1=1.7036,x2=-0.7961
%x0=[-1.5 1]; % x1=-1.7036 ,x2=0.7961


xx = fmincon(fun, x0, [],[],[],[],lb,ub,[],opts)

%% gardient vector and Hessian
x=sym('x',[1 2]);
G=simplify(gradient((4-2.1*x(1).^2+(x(1).^4)/3).*x(1).^2 + x(1).*x(2) + (-4+4*x(2).^2).*x(2).^2,x))
H=simplify(hessian((4-2.1*x(1).^2+(x(1).^4)/3).*x(1).^2 + x(1).*x(2) + (-4+4*x(2).^2).*x(2).^2,x));

% optimal solution x0=[-1 -1];
x1=-0.0899 ; x2=0.7126;
% gardient vector
eval(G)
% Hessian, eigenvalues>=0---> Hessian>=0
eval(H);eig(eval(H))

