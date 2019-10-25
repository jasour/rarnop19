% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 8: Nonlinear Robust Optimization
%% Robust Set X_R: Set of all design variable “x”  that satisfies “safety/design constraints” for all possible values of uncertainty “q”, i.e., g(x,q)<=0
% X_R = {x:  g(x,q)<=0 , for all values of uncertainty q } 
% Inner approximation X^d_R={x: p_d(x)<=0}
%
% Robust Control Set for Uncertain Nonlinear Dynamical System
% X^d_R={u(k): achieve control objectives for all possible values of uncertainties }
%%

clc;clear all;close all


% 2d: order of polynomial p_d(x) Robust Set = {x: p_d(x) <= 0}
d=4;d_sos=d/2;

% variables: design:x uncertainty : q
nx=1;nq=4;
x=sdpvar(1,nx);
q=sdpvar(1,nq);

% polynomial p_d(x) with unknown coefficients coef
vpow=[];for k = 0:2*d; vpow = [vpow;genpow(nx,k)]; end % powers of monomials 
coef=sdpvar(size(vpow,1),1); %coefficients 
p_d=coef'*(x.^vpow(:,1)); % polynomial p_d(x) 

% moments of lebesgue measure on [-1 1]:  yx(a)= \int x^a dx
% Integral of p_d(x) : replace x^a with yx(a)
yx=[2];for i=1:2*d ;yx(i+1,1)=(1/1)*(((1)^(i+1) - (-1)^(i+1))/(i+1));end 
Int_p_d=coef'*(yx(vpow(:,1)+1));

% nonlinear dynamics of system, 
% design variable: control input  u
% source of uncertainties at time k: (x1,x2,x3,w) states and disturbance at time k
 u=x(1); x1=q(1); x2=q(2); x3=q(3); w=q(4);
 f1=w*x2;
 f2=x1*x3;
 f3=1.2*x1-0.5*x2+x3+2*u;
 
% set K: uncertain control objective g(x,q) <=0
% states at time (k+1) \in neighbourhood of the given way-point (0,0,0.9)
g=-(1^2-(f1-0).^2/0.2^2-(f2-0).^2/0.2^2-(f3-0.9).^2/0.4^2)

% Set of all values of design variable u, i.e., u in X=[-1 1]
X=(u+1)*(1-u);

% Set of all values of uncertain parameter: states and disturbance at time k
Q=0.1^2-x1^2-x2^2-x3^2-w^2;

% sos polynomials
[s1,c1] = polynomial([x1;x2;x3;w],2*d_sos);


% SOS Conditions: p_d(x) > g(x,q) for all q in Q and x in X 
% In theory we need sos(p_d-g s1*Q -s2*X)
 F = [sos(p_d-g-s1*Q), sos(s1)];

% mosek solver
ops = sdpsettings('solver','mosek');

[sol,v,Q]=solvesos(F, Int_p_d,[],[c1;coef]);


%% Results

% obtained p_d(x)
syms x;
pp_d=value(coef)'*(x.^vpow(:,1));

%% Plots

figure
% p_d(x)
x=[-1:0.01:1];
plot(x,eval(pp_d),'LineWidth',5);hold on
plot(x,zeros(size(x)),'b','LineWidth',1)
xlabel('$x$','Interpreter','latex', 'FontSize',40);ylabel('$p_d(x)$','Interpreter','latex', 'FontSize',40);
str1 = '$\chi^d_{R}=\{x: p_d(x) \leq 0 \}$';text(0.4,0.1,str1,'Interpreter','latex','FontSize',30)
xlim([0 1]);ylim([-0.5 1]); grid on;hold on


%X^d_R= {x : p_d(x) <=0 }
x=[-1:0.002:1];ind=find(eval(pp_d)<=0); x=x(ind);xh=[x;eval(pp_d)]';
%plot([[xh(:,1)]'; [xh(:,1)]'],[ [xh(:,2)]';[0*ones(size(xh(:,1)))]'  ],'k--*','LineWidth',2 )
plot([[xh(:,1)]'],[ [0*ones(size(xh(:,1)))]'  ],'k--*','LineWidth',2 )

ax = gca; ax.FontSize = 30;
xlabel('$x$','Interpreter','latex', 'FontSize',40);ylabel('$p_d(x)$','Interpreter','latex', 'FontSize',40);
