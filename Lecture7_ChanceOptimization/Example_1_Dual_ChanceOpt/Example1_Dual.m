% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 7: Chance Constrained/Chance Optimization

%% This code obtains the Dual optimal values of the SDP in (3.7)
%  which approximates the chance optimization problem in (1.2)
%  e.q, Find x such that Probability of set { q : p_j(x,q)>0, j=1,...l } becomes maximum 
%       where q: random variable with probability distribution muq. 
%
%  "SEMIDEFINITE PROGRAMMING FOR CHANCE CONSTRAINED OPTIMIZATION OVER SEMIALGEBRAIC SETS"
% SIAM J. OPTIM. Vol. 25, No. 3, pp. 1411–1440
% ASHKAN. M. JASOUR, N. S. AYBAT, AND C. M. LAGOA
%% Dual of Moment SDP (SOS Optimization)

clc;clear;close all

% Relaxarion order:  polynomial order of 2d
d=10;d_sos=d;

% design and uncertain parameters
nx=1;nq=1;
x=sdpvar(1,nx);
q=sdpvar(1,nq);

% Upper bound of probability 
sdpvar Beta

% polynomial W(x,q) of order 2d
vpow=[];for k = 0:d; vpow = [vpow;genpow(2,k)]; end % monomials
coef=sdpvar(size(vpow,1),1); %coefficients 
W=coef'*(x.^vpow(:,1).*q(1).^vpow(:,2)); % polynomial W(x,q) 

% yq_i: given moments of measure muq_i, i=1,...,nq
% moments of uniform distribution on [-1,1]
yq_1=[1];for i=1:d ;yq_1(i+1,1)=(1/2)*(((1)^(i+1) - (-1)^(i+1))/(i+1));end 

% To calculate the integral of W(x,q) with respect to probability measure
% of q, in W(x,q) replace q^a with a-th moment of q
W_Int=coef'*(x.^vpow(:,1).*yq_1(vpow(:,2)+1));

% set K
K=0.5*q*(q^2+(x-0.5)^2)-(q^4+q^2*(x-0.5)^2+(x-0.5)^4); 

% set X
X=(x+1)*(1-x);

% sos polynomials
[s1,c1] = polynomial([x q],d_sos);
[s2,c2] = polynomial([x],d_sos);

% SOS constraints 
F = [sos(W-1-[s1]*K), sos(s1), sos(W), sos(Beta-W_Int-[s2]*X), sos(s2) ];

% SDP solver
ops = sdpsettings('solver','mosek');

% SOS program
[sol,v,Q]=solvesos(F, Beta,[],[c1;c2;Beta;coef]);

%% results
% obtained Beta 
BBeta=value(Beta);

x=sym('x',[1 nx]); q=sym('q',[1 nq]);
% obtained polynomial W(x,q)
WW=value(coef)'*(x(1).^vpow(:,1).*q(1).^vpow(:,2));
% obtained Integral of polynomial W(x,q) with respect to probability measure
WW_Int=value(coef)'*(x(1).^vpow(:,1).*yq_1(vpow(:,2)+1));

%% Plots

% plot W
[x1,q1]=meshgrid([-0.9:0.01:1],[-0.9:0.01:0.9]);
surf(x1,q1,eval(WW),'FaceColor','red','FaceAlpha',0.85,'EdgeColor','none','FaceLighting','phong');hold on;grid on;
surf(x1,q1,ones(size(eval(WW))),'FaceColor','blue','FaceAlpha',0.5,'EdgeColor','none','FaceLighting','phong')
axis tight;hidden off;ylabel('q');xlabel('x');hold on;xlim([-1,1]);ylim([-1,1]);zlim([0,3]);
camlight; lighting gouraud,hold on
title('$W(x,\omega)$','Interpreter','latex', 'FontSize',31);
zlim([0 2])

% plot set K1
for x=0:0.05:1; for q=0:0.05:1
        if (0.5*q*(q^2+(x-0.5)^2)-(q^4+q^2*(x-0.5)^2+(x-0.5)^4)) >=0; plot3([x x],[q q],[0 1],'k-*');end
end;end
xlabel('$x$','Interpreter','latex', 'FontSize',31);ylabel('$\omega$','Interpreter','latex', 'FontSize',31)
str1 = '$ \mathcal{K} $';text(0.94,0.2,str1,'HorizontalAlignment','right','Interpreter','latex','FontSize',30) 

% plot Integral of W 
figure
x1=[-0.9:0.01:0.9];
plot(x1,eval(WW_Int),x1,BBeta*ones(size(x1)),'LineWidth',5);grid;hold on
%title('$ \int \mathcal{P}^d_{\mathcal{W}}(x,q) d\mu_q$, \ \ \beta$','Interpreter','latex', 'FontSize',31);
xlabel('$x$','Interpreter','latex', 'FontSize',31);set(gca,'fontsize',20)
str2 = '$ x^* =0.5 $';text(0.5,0.3,str2,'HorizontalAlignment','right','Interpreter','latex','FontSize',30)         
str3 = '$ \mathbf{P_d^*} = \beta $';text(0.5,0.55,str3,'HorizontalAlignment','right','Interpreter','latex','FontSize',30)  
str4 = '$ \int {\mathcal{W}}(x,\omega) d\mu_{\omega}$';text(-0.3,0.4,str4,'HorizontalAlignment','right','Interpreter','latex','FontSize',30)  

ylim([0.2 0.6])