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
%%
clc;clear;close all

% Risk level delta
Delta=0.2;

% Relaxarion order:  h(x): polynomial order of 2d_h
d_h=8; % 20
% Relaxarion order:  W(x,q): polynomial order of 2d_w
d_w=8;
% Relaxarion order:  SOS polynommials 2d_sos
d_sos=5;

% design and uncertain parameters
nx=1;nq=1;
x=sdpvar(1,nx);
q=sdpvar(1,nq);

% polynomial W(x,q) of order 2d_w
vpow=[];for k = 0:2*d_w; vpow = [vpow;genpow(nx+nq,k)]; end % monomials
coef=sdpvar(size(vpow,1),1); %coefficients 
W=coef'*(x(1).^vpow(:,1).*q(1).^vpow(:,2)); % polynomial W(x,q) 

% polynomial h(x) of order 2d_h
h_vpow=[];for k = 0:2*d_h; h_vpow = [h_vpow;genpow(nx,k)]; end % monomials
h_coef=sdpvar(size(h_vpow,1),1); %coefficients 
h=h_coef'*(x(1).^h_vpow(:,1)); % polynomial h(x) 

% q has normal distribution with uncertain mean and sigma
% mean in [-0.1, 0.1]; sigma in [0.1 0.3]
% moments of q are polynomial function of mean and sigma
%Polynomial moments of uncertain normal distribution;
mom_poly;yq_1=py;

% Uncertain parameters of Normal distribution
% m: mean in A1=[-0.1, 0.1],  % ss: sigma in A2=[0.1, 0.3]
A1=0.1^2-mm^2;
A2=(0.3-ss)*(ss-0.1);

X=1-x(1)^2;

% moments of lebesgue measure on [-1,1] : measure with density function 1
yx_1=[2];for i=1:2*d_h ;yx_1(i+1,1)=(1/1)*(((1)^(i+1) - (-1)^(i+1))/(i+1));end 

% Int W(x,q) pr(q)dq
% To calculate the integral of W(x,q) with respect to probability measure
% of q, in W(x,q) replace q^a with a-th moment of q. 
W_Int=coef'*(x(1).^vpow(:,1).*yq_1(vpow(:,2)+1));

% Int h(x) dx
% To calculate the integral of W(x,q) with respect to lebesgue measure (dx)
%, in h(x) replace x^a with a-th moment of lebesgue measure
h_Int=h_coef'*(yx_1(h_vpow(:,1)+1));

% Safety Constraint, set K={(x,q):  g(x,q)>0}
% K_bar: complement set K_bar={(x,q):  g(x,q)<=0}
 K_bar=-(0.8^2-x(1)^2-q(1)^2); 

% sos polynomials
[s1,c1] = polynomial([mm ss],2*d_sos);
[s2,c2] = polynomial([mm ss],2*d_sos);
[s3,c3] = polynomial([x q],2*d_sos);

% SOS constraints 
% h(x) >= int W(x,q)pr(q)dq  for all uncertain parameters and all x
% W(x,q) >= 1 on set K_bar
% W(x,q) >=0 
F = [sos(h-W_Int-s1*A1-s2*A2), sos(W-1-s3*K_bar), sos(W),sos(s1),sos(s2),sos(s3)];
% In theory we need SOS(h-W_Int-s1*A1-s2*A2-s*X) where [s,c] = polynomial([x],2*d_sos) and X=[-1,1];


% SDP solver
ops = sdpsettings('solver','mosek');

% SOS program, min int h(x)dx s.t SOS conditions
[sol,v,Q]=solvesos(F, h_Int,[],[c1;c2;c3;coef;h_coef]);

%% Results

x=sym('x',[1 nx]); q=sym('q',[1 nq]); %syms mm ss
% obtained polynomial W(x,q)
WW=value(coef)'*(x(1).^vpow(:,1).*q(1).^vpow(:,2));
% obtained polynomial h(x)
hh=value(h_coef)'*(x(1).^h_vpow(:,1));

%% Plots

% plot W
[x1,q1]=meshgrid([-1:0.01:1],[-1:0.01:1]);hold on
surf(x1,q1,eval(WW),'FaceColor','red','FaceAlpha',0.85,'EdgeColor','none','FaceLighting','phong');hold on;grid on;
surf(x1,q1,ones(size(eval(WW))),'FaceColor','blue','FaceAlpha',0.5,'EdgeColor','none','FaceLighting','phong')
axis tight;hidden off;ylabel('q');xlabel('x');hold on;xlim([-1,1]);ylim([-1,1]);zlim([0,3]);
camlight; lighting gouraud,hold on
title('$W(x,\omega)$','Interpreter','latex', 'FontSize',31);
zlim([0 2])

% plot set K_bar
for x=-1:0.1:1; for q=-1:0.1:1
        if (0.8^2-x^2-q^2) >=0; plot3([x x],[q q],[0 1],'k-*');end
end;end
xlabel('$x$','Interpreter','latex', 'FontSize',31);ylabel('$\omega$','Interpreter','latex', 'FontSize',31)
str1 = '$ \mathcal{K} $';text(0.94,0.2,str1,'HorizontalAlignment','right','Interpreter','latex','FontSize',30) 


% plot h(x) 
figure
x1=[-1:0.01:1];
plot(x1,eval(hh),'color','red','LineWidth',5);grid;hold on
plot(x1,Delta*ones(size(x1)),'--','color','red','LineWidth',5);grid;hold on
xlabel('$x$','Interpreter','latex', 'FontSize',31);set(gca,'fontsize',20)
str2 = '$ \Delta $';text(0.6,0.2,str2,'Interpreter','latex','FontSize',30,'color','red')         
str3 = '$ h(x)$';text(-0.6,0.4,str3,'Interpreter','latex','FontSize',30,'color','red')  
ylim([-0.99 0.99]);ylim([0 1.5]); grid on
