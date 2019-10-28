% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 8: Nonlinear Robust Optimizaion
%% Distributionally Robust Chance Constrained Set 
% uncertain parameter q: parametric probability distribution (probabilistic sense), parameters are in
% bounded set(robust sense)
% x: design variable
% g(x,q)>=0 : safety Constraints
% Distributionally Robust Chance Constrained Set:
% X_DR = {all x:   probability { g(x,q)>=0 }>=1-Delta, for all probability distribution of uncertainty parameter q } 
% Inner approximation X^d_DR={x: p_d(x)<=delta}
%% Safe Control of Uncertain Nonlinear dynamics, 
% x: control input u, 
% q: uncertainties x1(k), x2(k), x3(k), omega(k)
% q1~ Normal(m1,sigma1)
% q2~ Normal(m2,sigma2)
% q3~ Normal(m3,sigma3)
% q4~ Normal(m4,sigma4)
% (m1,m2,m3,m4) \in Am={ 0.1-m1^2-m2^2-m3^2-m4^2 >=0 }
% (sigma1,sigma2,sigma3,sigma4) \in As={ 0.1-(sigma1-0.2)^2-(sigma2-0.2)^2-(sigma3-0.2)^2-(sigma4-0.2)^2 >=0 }
%%
clc;clear;close all

% Risk level delta
Delta=0.2;

% Relaxarion order:  h(x): polynomial order of 2d_h
d_h=6; % 6
% Relaxarion order:  W(x,q): polynomial order of 2d_w
d_w=4; %4
% Relaxarion order:  SOS polynommials 2d_sos
d_sos=2;

% design and uncertain parameters
nx=1;nq=4;
x=sdpvar(1,nx);
q=sdpvar(1,nq);

% polynomial W(x,q) of order 2d_w
vpow=[];for k = 0:2*d_w; vpow = [vpow;genpow(nx+nq,k)]; end % monomials
coef=sdpvar(size(vpow,1),1); %coefficients 
W=coef'*(x(1).^vpow(:,1).*q(1).^vpow(:,2).*q(2).^vpow(:,3).*q(3).^vpow(:,4).*q(4).^vpow(:,5)); % polynomial W(x,q) 

% polynomial h(x) of order 2d_h
h_vpow=[];for k = 0:2*d_h; h_vpow = [h_vpow;genpow(nx,k)]; end % monomials
h_coef=sdpvar(size(h_vpow,1),1); %coefficients 
h=h_coef'*(x(1).^h_vpow(:,1)); % polynomial h(x) 

% normal distributions with uncertain mean and sigma
% mean in set Am1; sigma in set As2
% moments are polynomial function of mean and sigma
%Polynomial moments of uncertain normal distribution;
mom_poly1;yq_1=py;
mom_poly2;yq_2=py;
mom_poly3;yq_3=py;
mom_poly4;yq_4=py;

% Uncertain parameters of Normal distribution
% m: mean in Am1
Am=0.1^2-mm1^2-mm2^2-mm3^2-mm4^2; 
% sigma: in As1
As=0.1^2-(ss1-0.2)^2-(ss2-0.2)^2-(ss3-0.2)^2-(ss4-0.2)^2;


% moments of lebesgue measure on [-1,1] : measure with density function 1
yx_1=[2];for i=1:2*d_h ;yx_1(i+1,1)=(1/1)*(((1)^(i+1) - (-1)^(i+1))/(i+1));end 

% Int W(x,q) pr(q)dq
% To calculate the integral of W(x,q) with respect to probability measure
% of q, in W(x,q) replace q^a with a-th moment of q. 
W_Int=coef'*(x(1).^vpow(:,1).*yq_1(vpow(:,2)+1).*yq_2(vpow(:,3)+1).*yq_3(vpow(:,4)+1).*yq_4(vpow(:,5)+1));

% Int h(x) dx
% To calculate the integral of W(x,q) with respect to lebesgue measure (dx)
%, in h(x) replace x^a with a-th moment of lebesgue measure
h_Int=h_coef'*(yx_1(h_vpow(:,1)+1));

% Safety Constraint, set K={(x,q):  g(x,q)>0}
% K_bar: complement set K_bar={(x,q):  g(x,q)<=0}
% Uncertain Nonlinear dynamics
 u=x(1); x1=q(1); x2=q(2); x3=q(3); w=q(4);
f1=0.2*w*x2;
f2=x1*x3;
f3=1.2*x1-0.5*x2+x3+2*u;
K_bar=1^2-(f1-0.1).^2/0.2^2-(f2-0.1).^2/0.2^2-(f3-0.4).^2/0.3^2;


% sos polynomials
[Sm,Cm] = polynomial([mm1 mm2 mm3 mm4],2*d_sos);
[Ss,Cs] = polynomial([ss1 ss2 ss3 ss4],2*d_sos);


[s,c] = polynomial([x q],2*d_sos);

% SOS constraints 
% h(x) >= int W(x,q)pr(q)dq  for all uncertain parameters and all x
% W(x,q) >= 1 on set K_bar
% W(x,q) >=0 
F = [sos(h-W_Int-Sm*Am-Ss*As), sos(W-1-s*K_bar), sos(W),sos(Sm),sos(Ss),sos(s)];
% In theory we need SOS(h-W_Int-s1*A1-s2*A2-s*X) where [s,c] = polynomial([x],2*d_sos) and X=[-1,1];


% SDP solver
ops = sdpsettings('solver','mosek');

% SOS program, min int h(x)dx s.t SOS conditions
[sol,v,Q]=solvesos(F, h_Int,[],[Cm;Cs;c;coef;h_coef]);

%% Results

x=sym('x',[1 nx]); q=sym('q',[1 nq]); %syms mm ss
% obtained polynomial W(x,q)
%WW=value(coef)'*(x(1).^vpow(:,1).*q(1).^vpow(:,2));
% obtained polynomial h(x)
hh=value(h_coef)'*(x(1).^h_vpow(:,1));

%% Plots

% plot h(x) 
figure
x1=[-1:0.01:1];
plot(x1,eval(hh),'color','red','LineWidth',5);grid;hold on
plot(x1,Delta*ones(size(x1)),'--','color','red','LineWidth',5);grid;hold on
xlabel('$u$','Interpreter','latex', 'FontSize',31);set(gca,'fontsize',20)
str2 = '$ \Delta $';text(0.6,0.2,str2,'Interpreter','latex','FontSize',30,'color','red')         
str3 = '$ h(u)$';text(-0.6,0.4,str3,'Interpreter','latex','FontSize',30,'color','red')  
ylim([-0.99 0.99]);ylim([0 1.5]); grid on
