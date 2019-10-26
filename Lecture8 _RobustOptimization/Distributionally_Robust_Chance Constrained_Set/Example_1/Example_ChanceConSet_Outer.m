% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 7: Chance Constrained/Chance Optimization

%% Chance Constrained Set (SOS Optimization)

clc;clear;close all

% Risk level delta
Delta=0.2;
% Relaxarion order:  polynomial order of 2d
d=8;d_sos=d;

% design and uncertain parameters
nx=1;nq=1;
x=sdpvar(1,nx);
q=sdpvar(1,nq);

% polynomial W(x,q) of order 2d
vpow=[];for k = 0:2*d; vpow = [vpow;genpow(nx+nq,k)]; end % monomials
coef=sdpvar(size(vpow,1),1); %coefficients 
W=coef'*(x(1).^vpow(:,1).*q(1).^vpow(:,2)); % polynomial W(x,q) 

% yq_i: given moments of measure muq_i, i=1,...,nq
% moments of Norma distribution mean=mm, sigma=ss
mm=0;ss=0.2;
yq_1=[1]; for n=1:2*d; yq_1=[yq_1;ss^n*(-j*sqrt(2))^n*kummerU(-n/2, 1/2, -1/2*(mm/ss)^2)];    end



% moments of lebesgue measure on [-1,1] : measure with density function 1
yx_1=[2];for i=1:2*d ;yx_1(i+1,1)=(1/1)*(((1)^(i+1) - (-1)^(i+1))/(i+1));end 

% To calculate the integral of W(x,q) with respect to probability measure
% of q, in W(x,q) replace q^a with a-th moment of q and x^a with a-th moment of lebesgue measure
W_Int=coef'*(yx_1(vpow(:,1)+1).*yq_1(vpow(:,2)+1));

% set K
 K=0.8^2-x(1)^2-q(1)^2; 


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
WW=value(coef)'*(x(1).^vpow(:,1).*q(1).^vpow(:,2));
% obtained Integral of polynomial W(x,q) with respect to probability measure
WW_Int=value(coef)'*(x(1).^vpow(:,1).*yq_1(vpow(:,2)+1));

%% Plots

% plot W
[x1,q1]=meshgrid([-1:0.01:1],[-1:0.01:1]);hold on
surf(x1,q1,eval(WW),'FaceColor','red','FaceAlpha',0.85,'EdgeColor','none','FaceLighting','phong');hold on;grid on;
surf(x1,q1,ones(size(eval(WW))),'FaceColor','blue','FaceAlpha',0.5,'EdgeColor','none','FaceLighting','phong')
axis tight;hidden off;ylabel('q');xlabel('x');hold on;xlim([-1,1]);ylim([-1,1]);zlim([0,3]);
camlight; lighting gouraud,hold on
title('$W(x,\omega)$','Interpreter','latex', 'FontSize',31);
zlim([0 2])

% plot set K1
for x=-1:0.1:1; for q=-1:0.1:1
        if (0.8-x^2-q^2) >=0; plot3([x x],[q q],[0 1],'k-*');end
end;end
xlabel('$x$','Interpreter','latex', 'FontSize',31);ylabel('$\omega$','Interpreter','latex', 'FontSize',31)
str1 = '$ \mathcal{K} $';text(0.94,0.2,str1,'HorizontalAlignment','right','Interpreter','latex','FontSize',30) 

% plot Integral of W 
figure
x1=[-1:0.01:1];
plot(x1,eval(WW_Int),'color','red','LineWidth',5);grid;hold on
%plot([-1:0.01:1],1-Delta*ones(size(x1)),'--','color','red','LineWidth',3)
xlabel('$x$','Interpreter','latex', 'FontSize',31);set(gca,'fontsize',20)
str3 = '$ \int {\mathcal{W}}(x,\omega) d\mu_{\omega}$';text(-0.3,1.1,str3,'HorizontalAlignment','right','Interpreter','latex','FontSize',30,'color','red')  
%str3 = '$ \Delta$';text(-0.3,0.2,str3,'HorizontalAlignment','right','Interpreter','latex','FontSize',30,'color','red')  

ylim([0 1.5])
