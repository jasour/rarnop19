% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 3: Sum Of Squares Based SDP For Nonlinear Optimization
%% Moment Based SDP for Constrained Nonlinear Optimization : min_x p(x) s.t g_i(x)>=0  , i=1,...,m

clc;clear all;close all

% SDP solvers
mset clear; warning('off','YALMIP:strict')
mset('yalmip',true);mset(sdpsettings('solver','mosek')); % SDP sovers: mosek, sedumi, sdpt3,...

% variables x1 x2
mpol x1 x2
% objective function p(x1,x)
p = -x1;
% Constraints g_i(x)>=0, i=1,2,3
g=[3-2*x2-x1^2-x2^2;-x1-x2-x1*x2;1+x1*x2];

% d: relaxation order. SDP will be based on the moments up to order 2d
d=1;

% Generate moment SDP of order 2d 
P = msdp(min(p),g>=0,d);

% Solve Moment SDP
[status,obj] = msol(P)
 
 %% Results
 %% status==1: the moment SDP is solved successfully and Rank conditions are satisfied. Hence, GloptiPoly can extract the global optimal solutions.
 %% status==0: the moment SDP is solved successfully But Rank conditions are Not satisfied. Hence, GloptiPoly can NOT extract the global optimal solutions. Increase the relaxation order d.
 %% status==-1:  : moment SDP could NOT be solved (unbounded SDP).


if status==1
    
% extracted global optimal solution    
xx1=double(x1)
xx2=double(x2)

end



 