% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 4: Measure and Moments Based SDP For Nonlinear Optimization
%% Moment Based SDP for Unconstrained Nonlinear Optimization : min_x p(x)

clc;clear all;close all

% SDP solvers
mset clear; warning('off','YALMIP:strict')
mset('yalmip',true);mset(sdpsettings('solver','mosek')); % SDP sovers: mosek, sedumi, sdpt3,...

% variables x1 x2
mpol x1 x2
% objective function p(x1,x)
p=3+2*x1+2*x2+3*x1^2+2*x1*x2+3*x2^2+x1^4+x2^4;

% Generate moment SDP to miniize p(x)
P = msdp(min(p));

% Solve Moment SDP
[status,obj] = msol(P)
 
 %% Results
 %% status==1: the moment SDP is solved successfully and Rank conditions are satisfied. Hence, GloptiPoly can extract the global optimal solutions.
 %% status==0: the moment SDP is solved successfully But Rank conditions are Not satisfied. Hence, GloptiPoly can NOT extract the global optimal solutions. 
 %% status==-1:  : moment SDP could NOT be solved (unbounded SDP).


if status==1
    
% extracted global optimal solution    
xx1=double(x1)
xx2=double(x2)

end



 