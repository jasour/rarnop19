% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 7: Chance Constrained/Chance Optimization

%% Chance Constrained Set (SOS Optimization)

clc;clear;close all

% Risk level delta
Delta=0.2;
% Relaxarion order:  polynomial order of 2d
d=3;d_sos=d;

% design and uncertain parameters
nx=1;nq=4;
x=sdpvar(1,nx);
q=sdpvar(1,nq);

% polynomial W(x,q) of order 2d
vpow=[];for k = 0:2*d; vpow = [vpow;genpow(nx+nq,k)]; end % monomials
coef=sdpvar(size(vpow,1),1); %coefficients 
%W=coef'*(x.^vpow(:,1).*q(1).^vpow(:,2)); % polynomial W(x,q) 
W=coef'*(x.^vpow(:,1).*q(1).^vpow(:,2).*q(2).^vpow(:,3).*q(3).^vpow(:,4).*q(4).^vpow(:,5));

% yq_i: given moments of measure muq_i, i=1,...,nq
% moments of uniform distribution on [-1,1]
l1=-0.1;u1=0.1;yq_1=[1];for i=1:2*d ;yq_1(i+1,1)=(1/(u1-l1))*(((u1)^(i+1) - (l1)^(i+1))/(i+1));end 
l2=-0.1;u2=0.1;yq_2=[1];for i=1:2*d ;yq_2(i+1,1)=(1/(u2-l2))*(((u2)^(i+1) - (l2)^(i+1))/(i+1));end 
l3=-0.1;u3=0.1;yq_3=[1];for i=1:2*d ;yq_3(i+1,1)=(1/(u3-l3))*(((u3)^(i+1) - (l3)^(i+1))/(i+1));end 
% moments of Beta(alpha,beta) distribution 
alpha=2;beta=5; yq_4=[1];for k=1:2*d; yq_4=[yq_4,(alpha+k-1)/(alpha+beta+k-1)*yq_4(end) ]; end;  yq_4=yq_4';

% moments of lebesgue measure on [-1,1] : measure with density function 1
yx_1=[2];for i=1:2*d ;yx_1(i+1,1)=(1/1)*(((1)^(i+1) - (-1)^(i+1))/(i+1));end 

% To calculate the integral of W(x,q) with respect to probability measure
% of q, in W(x,q) replace q^a with a-th moment of q and x^a with a-th moment of lebesgue measure
%W_Int=coef'*(yx_1(vpow(:,1)+1).*yq_1(vpow(:,2)+1));
W_Int=coef'*(yx_1(vpow(:,1)+1).*yq_1(vpow(:,2)+1).*yq_2(vpow(:,3)+1).*yq_3(vpow(:,4)+1).*yq_4(vpow(:,5)+1));

% set K
u=x(1); x1=q(1); x2=q(2); x3=q(3); w=q(4);
 f1=0.2*w*x2;
 f2=x1*x3;
 f3=1.2*x1-0.5*x2+x3+2*u;
 K=1^2-(f1-0).^2/0.03^2-(f2-0).^2/0.02^2-(f3-1).^2/0.4^2;

% sos polynomials
[s1,c1] = polynomial([x q],2*d_sos);

% SOS constraints 
F = [sos(W-1-[s1]*K), sos(s1), sos(W) ];

% SDP solver
ops = sdpsettings('solver','mosek');

% SOS program
[sol,v,Q]=solvesos(F, W_Int,[],[c1;coef]);

%% results

x=sym('x',[1 nx]); q=sym('q',[1 nq]);
% obtained polynomial W(x,q)
%WW=value(coef)'*(x(1).^vpow(:,1).*q(1).^vpow(:,2));
WW=value(coef)'*(x(1).^vpow(:,1).*q(1).^vpow(:,2).*q(2).^vpow(:,3).*q(3).^vpow(:,4).*q(4).^vpow(:,5));

% obtained Integral of polynomial W(x,q) with respect to probability measure
%WW_Int=value(coef)'*(x(1).^vpow(:,1).*yq_1(vpow(:,2)+1));
WW_Int=value(coef)'*(x(1).^vpow(:,1).*yq_1(vpow(:,2)+1).*yq_2(vpow(:,3)+1).*yq_3(vpow(:,4)+1).*yq_4(vpow(:,5)+1));

%% Plots

% plot Integral of W 
figure
x1=[-0.9:0.01:1];
plot(x1,eval(WW_Int),x1,1-Delta*ones(size(x1)),'LineWidth',5);grid;hold on
%title('$ \int \mathcal{P}^d_{\mathcal{W}}(x,q) d\mu_q$, \ \ \beta$','Interpreter','latex', 'FontSize',31);
xlabel('$u$','Interpreter','latex', 'FontSize',31);set(gca,'fontsize',20)
str2 = '$ 1-\Delta $';text(0.1,0.9,str2,'HorizontalAlignment','right','Interpreter','latex','FontSize',30)         
str3 = '$ \int {\mathcal{W}}(x,\omega) d\mu_{\omega}$';text(-0.3,0.1,str3,'HorizontalAlignment','right','Interpreter','latex','FontSize',30)  
ylim([0 1.5])

pause(0.1)

% Monte Carlo Probability Curve
disp('Generating Monte Carlo Probability Curve')
Example_3_MonteCarlo;
str4 = 'Monte Carlo Probability Curve';text(0.8,0.3,str4,'HorizontalAlignment','right','Interpreter','latex','FontSize',30)  
ylim([0 1.5])