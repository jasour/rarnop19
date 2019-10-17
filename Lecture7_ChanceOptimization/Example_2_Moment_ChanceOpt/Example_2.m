% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 7: Chance Constrained/Chance Optimization

%% This code obtains the optimal values of the SDP in (3.7)
%  which approximates the chance optimization problem in (1.2)
%  e.q, Find x such that Probability of set { q : p_j(x,q)>0, j=1,...l } becomes maximum 
%       where q: random variable with probability distribution muq. 
%
%  "SEMIDEFINITE PROGRAMMING FOR CHANCE CONSTRAINED OPTIMIZATION OVER SEMIALGEBRAIC SETS"
% SIAM J. OPTIM. Vol. 25, No. 3, pp. 1411–1440
% ASHKAN. M. JASOUR, N. S. AYBAT, AND C. M. LAGOA
%%

clc;clear;clear all; 

%% Problem Parameters
d =5; % relaxation order, looks for the moments up to order 2d
nx =1;% number of decision variables : x
nq= 1;% number of uncertain parameters : q

%% Mesures mu, mu_x, mu_s  mu_q     :   mu_s= mu_x*mu_q - mu
% mu defined in space (x,q)
% mu_x defined in space (x) (denoted by x2)
% mu_s defined in space (x,q) (denoted by (xs,qs))
% To distinguish measures, we use different variables to represent the related space.

% mu: measure supported on p_j>=0 j=1,...n_p, y: moments of mu,  
mpol('x',nx); mpol('q',nq); mu = meas([x;q]); y=mom(mmon([x;q],2*d)); 
% success set K ={ p_j(x,q)>0, j=1,...n_p }
K=-[((3*x(1))/2 + 1).^2/4 - ((3*x(1))/2 + 1).^3/4 + ((3*x(1))/2 + 1).^4/16 + (9*q(1).^2)/100 - 29/400];

% mu_x: measure, yx: moments of mu_x, 
mpol('x2',nx);mux= meas([x2]); yx=mom(mmon([x2],2*d)); 

% mu_s: slack measure, mu_s = mux*muq - mu, y_s: moments of mu s, 
mpol('x_s',nx); mpol('q_s',nq); mu_s = meas([x_s;q_s]); y_s=mom(mmon([x_s;q_s],2*d));

% yq_i: given moments of measure muq_i, i=1,...,nq
% moments of uniform distribution on [-1,1]
yq_1=[1];for i=1:2*d ;yq_1(i+1,1)=(1/2)*(((1)^(i+1) - (-1)^(i+1))/(i+1));end 
% yq=[yq_1,yq_2,...,yq_nq]
yq=[];for i=1:nq; yq=[yq,yq_1]; end

% yxq: moments of cross product measure mux*muq
yxq=  mom_xq(nx,nq,d,yq,yx);

%% Optimization
mset('yalmip',true);mset(sdpsettings('solver','mosek')); % Calls mosek SDP solver through Yalmip   

% Chance Optimization in measures
Opt=msdp(max(mass(mu)),mass(mux)==1,K>=0,y_s==(yxq - y),-1<=yx,yx<=1,d);

% Solves the moment SDP
msol(Opt);

%% Results
% moments of measure mu
y=double(mvec(mu)); 
% moment matrix of measure mu
M=double(mmat(mu)); 
% moments of measure mux
yx=double(mvec(mux)); 
% moment matrix of measure mux
Mx=double(mmat(mux)); 
%%
%% upper bound of obtained probability
Probailty =y(1); 

%%
%% Rank Test of Moment matrix of design variable x: If Rank(Md)=Rank(Md-1)= r : r Dirac measure : r global optimal solution. (Md is flat extension of Md-1)  

% Vector of monomials up to order d
d2= d;
B_d2=mmon([x2],d2);
% Moment Matrix of order d: Md=E[Bd*Bd']
M_d2=double(mom(B_d2*B_d2'));
% Rank of Md: nonzero eigenvalues
R_d2=eig(M_d2);

% Vector of monomials up to order d-1
d1= d2-1;
B_d1=mmon([x2],d1);
% Moment Matrix of order d: Md=E[Bd*Bd']
M_d1=double(mom(B_d1*B_d1'));
% Rank of Md: nonzero eigenvalues
R_d1=eig(M_d1);

%%
%% Extract x if rank condition is satisfied
x = extract_mom(yx,nx,d);

%%
 disp('Probailty:');Probailty
 disp('Decision parameter x:');x


