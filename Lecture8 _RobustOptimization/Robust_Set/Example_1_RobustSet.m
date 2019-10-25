% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 8: Nonlinear Robust Optimization
%% Robust Set X_R: Set of all design variable “x”  that satisfies “safety/design constraints” for all possible values of uncertainty “q”, i.e., g(x,q)<=0
% X_R = {x:  g(x,q)<=0 , for all values of uncertainty q } 
% Inner approximation X^d_R={x: p_d(x)<=0}

%%
clc;clear all;close all


% 2d: order of polynomial p_d(x) Robust Set = {x: p_d(x) <= 0}
d=1;d_sos=d;
% variables: design:x uncertainty : q
sdpvar x q 

% polynomial p_d(x) with unknown coefficients coef
vpow=[];for k = 0:2*d; vpow = [vpow;genpow(1,k)]; end % powers of monomials 
coef=sdpvar(size(vpow,1),1); %coefficients 
p_d=coef'*(x.^vpow(:,1)); % polynomial p_d(x) 

% moments of lebesgue measure on [-1 1]:  yx(a)= \int x^a dx
% Integral of p_d(x) : replace x^a with yx(a)
yx=[2];for i=1:2*d ;yx(i+1,1)=(1/1)*(((1)^(i+1) - (-1)^(i+1))/(i+1));end 
Int_p_d=coef'*(yx(vpow(:,1)+1));

% set K: uncertain safety constraint g(x,q) <=0
g=x^2+q^2-1; 

% Set of all values of design variable x, i.e., x in X=[-1 1]
X=(x+1)*(1-x);
% Set of all values of uncertain parameter q, i.e., q in Q=[-0.5 0.5]
Q=(q+0.5)*(0.5-q);

% sos polynomials
[s1,c1] = polynomial([x,q],2*d_sos);
[s2,c2] = polynomial([x,q],2*d_sos);

% SOS Conditions: p_d(x) > g(x,q) for all q in Q and x in X 
F = [sos(p_d-g-s1*X-s2*Q), sos(s1),sos(s2)];

% mosek solver
ops = sdpsettings('solver','mosek');
[sol,v,Q]=solvesos(F, Int_p_d,[],[c1;c2;coef]);

%% Results

% obtained p_d(x)
syms x;
pp_d=value(coef)'*(x.^vpow(:,1));

%% Plots

figure
% p_d(x)
x=[-1:0.01:1];
plot(x,eval(pp_d),x,zeros(size(x)),'LineWidth',5)
xlabel('$x$','Interpreter','latex', 'FontSize',40);ylabel('$p_d(x)$','Interpreter','latex', 'FontSize',40);
str1 = '$\chi^d_{R}=\{x: p_d(x) \leq 0 \}$';text(-0.9,0.1,str1,'Interpreter','latex','FontSize',30)

figure; hold on

% g(x,q)
[x,q]=meshgrid([-2:0.01:2],[-0.5:0.01:0.5]);
surf(x,q,x.^2+q.^2-1,'FaceColor','red','FaceAlpha',0.9,'EdgeColor','none','FaceLighting','phong');
[x,q]=meshgrid([-2:0.01:2],[-2:0.01:2]);
h=surf(x,q,zeros(size(x)),'FaceColor','blue','FaceAlpha',0.7,'EdgeColor','none','FaceLighting','phong');
hold on;grid on;camlight; lighting gouraud,hold on
str1 = '$g(x,\omega)$';text(0,0,1,str1,'Interpreter','latex','FontSize',30,'color','red')

% g(x,q*), q*: worst case uncertainty  q*=argmin max g(x,q) s.t. q in Q
x=[-2:0.01:2];
plot3(x,0.5*ones(size(x)),x.^2+0.5.^2-1,'--k','LineWidth',5); hold on
plot3(x,-0.5*ones(size(x)),x.^2+0.5.^2-1,'--k','LineWidth',5)
str1 = '$g(x,\omega^*)$';text(-1,1.2,-0.5,str1,'Interpreter','latex','FontSize',30,'color','black')

% p_d(x)
x=[-2:0.1:2];
plot3(x,0.5*ones(size(x)),eval(pp_d),'--y','LineWidth',8)
str1 = '$p_d(x)$';text(-1,1.5,1,str1,'Interpreter','latex','FontSize',30,'color','yellow')

%X^d_R= {x : p_d(x) <=0 }
x=[-1:0.1:1];ind=find(eval(pp_d)<=0); x=x(ind);xh=[x;eval(pp_d)]';
plot3([[xh(:,1)]'; [xh(:,1)]'],[[0.5*ones(size(xh(:,1)))]'; [0.5*ones(size(xh(:,1)))]'],[ [xh(:,2)]';[0*ones(size(xh(:,1)))]'  ],'k--*','LineWidth',2 )

ax = gca; ax.FontSize = 30;
zlim([-1 2]); xlabel('$x$','Interpreter','latex', 'FontSize',40);ylabel('$\omega$','Interpreter','latex', 'FontSize',40); 
zlabel('$g(x,\omega)$','Interpreter','latex', 'FontSize',40);
