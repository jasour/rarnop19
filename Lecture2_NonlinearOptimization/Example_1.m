% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 2: Nonlinear Optimization

%% 

clc;clear;close all;

opts = optimoptions(@fmincon,'Algorithm','interior-point','PlotFcn',{@optimplotx,@optimplotfval},'Display', 'iter');


fun = @(x)1-x^2+x^3+1/4*x^4-3/5*x^5+1/6*x^6;
%lower bound
lb=[-3];
%uper bound
ub=[3];
% initial point 
x0=3;
% optimal solution
x= fmincon(fun, x0, [],[],[],[],lb,ub,[],opts)

pause

% initial point 
x0=-3;
% optimal solution
x= fmincon(fun, x0, [],[],[],[],lb,ub,[],opts)