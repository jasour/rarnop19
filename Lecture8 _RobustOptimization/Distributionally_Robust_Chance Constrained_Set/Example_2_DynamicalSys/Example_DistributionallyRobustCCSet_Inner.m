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
% q1,q2,q3 ~ Uniform [-0.1, 0.1]^3
% q4 ~ Normal(m,sigma),  m \in [-0.1, 0.1],  sigma \in [0.1, 0.3]

clc;clear;close all

% Risk level delta
Delta=0.2;

% Relaxarion order:  h(x): polynomial order of 2d_h
d_h=4; % 20
% Relaxarion order:  W(x,q): polynomial order of 2d_w
d_w=4;
% Relaxarion order:  SOS polynommials 2d_sos
d_sos=4;
d=d_w;
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



% moments of uniform distribution on [-1,1]
l1=-0.1;u1=0.1;yq_1=[1];for i=1:2*d ;yq_1(i+1,1)=(1/(u1-l1))*(((u1)^(i+1) - (l1)^(i+1))/(i+1));end 
l2=-0.1;u2=0.1;yq_2=[1];for i=1:2*d ;yq_2(i+1,1)=(1/(u2-l2))*(((u2)^(i+1) - (l2)^(i+1))/(i+1));end 
l3=-0.1;u3=0.1;yq_3=[1];for i=1:2*d ;yq_3(i+1,1)=(1/(u3-l3))*(((u3)^(i+1) - (l3)^(i+1))/(i+1));end 
% normal distribution with uncertain mean and sigma
% mean in [-0.1, 0.1]; sigma in [0.1 0.3]
% moments are polynomial function of mean and sigma
%Polynomial moments of uncertain normal distribution;
mom_poly4;yq_4=py;

% Uncertain parameters of Normal distribution
% m: mean in A1=[-0.1, 0.1],  % ss: sigma in A2=[0.1, 0.3]
Am4=0.1^2-mm4^2; As4=(0.3-ss4)*(ss4-0.1);



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
[Sm4,Cm4] = polynomial([mm4 ss4],2*d_sos);
[Ss4,Cs4] = polynomial([mm4 ss4],2*d_sos);

[s,c] = polynomial([x q],2*d_sos);

% SOS constraints 
% h(x) >= int W(x,q)pr(q)dq  for all uncertain parameters and all x
% W(x,q) >= 1 on set K_bar
% W(x,q) >=0 
F = [sos(h-W_Int-Sm4*Am4-Ss4*As4), sos(W-1-s*K_bar), sos(W),sos(Sm4),sos(Ss4),sos(s)];
% In theory we need SOS(h-W_Int-s1*A1-s2*A2-s*X) where [s,c] = polynomial([x],2*d_sos) and X=[-1,1];


% SDP solver
ops = sdpsettings('solver','mosek');

% SOS program, min int h(x)dx s.t SOS conditions
[sol,v,Q]=solvesos(F, h_Int,[],[Cm4;Cs4;c;coef;h_coef]);

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
str2 = '$ \Delta $';text(0.6,0.3,str2,'Interpreter','latex','FontSize',30,'color','red')         
str3 = '$ h(u)$';text(-0.6,0.4,str3,'Interpreter','latex','FontSize',30,'color','red')  
ylim([-0.99 0.99]);ylim([0 1.5]); grid on
