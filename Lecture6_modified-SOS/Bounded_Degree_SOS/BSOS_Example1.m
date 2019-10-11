% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 3: Sum Of Squares Based SDP For Nonlinear Optimization
%% Bounded degree SOS based SDP for Constrained optimization : min_x p(x)

clear all;clc;

%  F(x)= x1^2    x2^2        x3^2 
pop.F = [2 0 0 1; 0 2 0 1; 0 0 2 1]; %[term1: pow(x1) pow(x2) pow(x3) coeff;]
pop.n = 3; % number of variables
pop.I = {1:3}; % variables in objective function

% Constraints
pop.G{1} = [1 0 0 1];% x1^2>=0
pop.G{2} = [0 1 0 1];% x2^2>=0
pop.G{3} = [0 0 1 1];% x3^2>=0
pop.J = {[1],[2],[3]}; %variables in Cons


pop.k=1; %  
pop.d=2; % 


sdp = gendata2(pop,'BSOS');
sol = csol(sdp,'sqlp');
psol = postproc(pop,sdp,sol);

clc
% optimal x
psol.obj
% rank of moment mat
psol.rnk
% totla cpu time
psol.ttot
