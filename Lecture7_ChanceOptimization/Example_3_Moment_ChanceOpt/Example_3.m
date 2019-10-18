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
d =2; % relaxation order, looks for the moments up to order 2d
nx =1;% number of decision variables : x
nq= 4;% number of uncertain parameters : q

%% Mesures mu, mu_x, mu_s  mu_q     :   mu_s= mu_x*mu_q - mu
% mu defined in space (x,q)
% mu_x defined in space (x) (denoted by x2)
% mu_s defined in space (x,q) (denoted by (xs,qs))
% To distinguish measures, we use different variables to represent the related space.

% mu: measure supported on p_j>=0 j=1,...n_p, y: moments of mu,  
mpol('x',nx); mpol('q',nq); mu = meas([x;q]); y=mom(mmon([x;q],2*d)); 

% success set K ={ p_j(x,q)>0, j=1,...n_p }
u=x(1); x1=q(1); x2=q(2); x3=q(3); w=q(4);
 f1=0.2*w*x2;
 f2=x1*x3;
 f3=1.2*x1-0.5*x2+x3+2*u;
 K=1^2-(f1-0).^2/0.02^2-(f2-0).^2/0.03^2-(f3-1).^2/0.4^2;

% mu_x: measure, yx: moments of mu_x, 
mpol('x2',nx);mux= meas([x2]); yx=mom(mmon([x2],2*d)); 

% mu_s: slack measure, mu_s = mux*muq - mu, y_s: moments of mu s, 
mpol('x_s',nx); mpol('q_s',nq); mu_s = meas([x_s;q_s]); y_s=mom(mmon([x_s;q_s],2*d));

% yq_i: given moments of measure muq_i, i=1,...,nq
% moments of uniform distribution on [-1,1]
l1=-0.1;u1=0.1;yq_1=[1];for i=1:2*d ;yq_1(i+1,1)=(1/(u1-l1))*(((u1)^(i+1) - (l1)^(i+1))/(i+1));end 
l2=-0.1;u2=0.1;yq_2=[1];for i=1:2*d ;yq_2(i+1,1)=(1/(u2-l2))*(((u2)^(i+1) - (l2)^(i+1))/(i+1));end 
l3=-0.1;u3=0.1;yq_3=[1];for i=1:2*d ;yq_3(i+1,1)=(1/(u3-l3))*(((u3)^(i+1) - (l3)^(i+1))/(i+1));end 
% moments of Beta(alpha,beta)
alpha=2;beta=5;yq_4=[1];for k=1:2*d; yq_4=[yq_4,(alpha+k-1)/(alpha+beta+k-1)*yq_4(end) ]; end; yq_4=yq_4';

% all moments   
yq=[yq_1,yq_2,yq_3,yq_4];

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
Probailty =y(1); disp('Probailty:');Probailty

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
%% The rank of momen matrix is approximately 1; Hence x=E[x]=yx(2)
% obtained decision variables x
Decision_x = yx(2:1+nx); 
disp('Decision parameter x:');Decision_x

