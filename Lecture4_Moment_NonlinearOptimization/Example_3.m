% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 4: Measure and Moments Based SDP For Nonlinear Optimization
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

% d: relaxation order. SDP will be based on the moments up to order 2d.
% 2d>= max ( deg(p),deg(g_i) )
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


%%
%% Results based on the moments
%%

if status ~=-1

% obtained moment vector 
y=double(mvec(meas)); 
% Moment matrix of the moments
M=double(mmat(meas)); 

%% Rank Test: If Rank(Md)=Rank(Md-1)= r : r Dirac measure : r global optimal solution. (Md is flat extension of Md-1)  

% Vector of monomials up to order d;
d2= d;
B_d2=mmon([x1 x2],d2);
% Moment Matrix of order d: Md=E[Bd*Bd']
M_d2=double(mom(B_d2*B_d2'));
% Rank of Md: nonzero eigenvalues
R_d2=eig(M_d2);

% Vector of monomials up to order d-max(dg);
d1= d-ceil(deg(g)/2);
B_d1=mmon([x1 x2],d1);
% Moment Matrix of order d: Md=E[Bd*Bd']
M_d1=double(mom(B_d1*B_d1'));
% Rank of Md: nonzero eigenvalues
R_d1=eig(M_d1);

end% status check

 