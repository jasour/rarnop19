% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 10: Probabilistic Safety Verification
%% Moment Based SDP for upper bound risk estimation: probability that random variable x satisfies nonlinear safety constraints
% max {int d(mu_bar)} s.t. mu_bar <= mu,  supp(mu_bar) \in {x: g(x)>=0}
% mu: given probability of x
% g(x)>=0 : safety constraint
% mu_bar: unknow measure,  int d(mu_bar) = risk
%%
clc;clear all;close all
% SDP solvers
mset clear; warning('off','YALMIP:strict')
mset('yalmip',true);mset(sdpsettings('solver','mosek')); % SDP sovers: mosek, sedumi, sdpt3,...


% relaxaition order
d = 25; 

% random variable x 
mpol x 
m = meas(x); 
% given moments of the x: uniform distribution [-1,1]
yx=[1];for i=1:d ;yx(i+1,1)=(1/2)*(((1)^(i+1) - (-1)^(i+1))/(i+1));end 

% safety constraint
g = x*(1/2-x);


% unknown moments of mu_bar 
y_bar = mom(mmon(x,d));

% slack variable xs
mpol xs
ms = meas(xs);
ys = mom(mmon(xs,d));

% moment SDP
P = msdp(max(mass(m)), g>=0, ys==yx-y_bar); 
% solve moment SDP relaxation
msol(P); 

% Result
y_bar = double(mvec(m)); 
% upper bound of risk
Probability=y_bar(1)
