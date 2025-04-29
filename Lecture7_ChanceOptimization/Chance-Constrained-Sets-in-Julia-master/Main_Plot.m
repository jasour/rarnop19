clc;clear;close all

Co = dlmread('Sol_Co.txt');
d=10; Delta=0.85;
yw1=[1];for i=1:2*d ;yw1(i+1,1)=(1/2)*(((1)^(i+1) - (-1)^(i+1))/(i+1));end 
vpow=[];for k = 2*d:-1:0; vpow = [vpow;genpow(2,k)]; end % monomials
%% Plots
syms x1 w1

% obtained polynomial P(x1,w1)
P=Co'*(x1.^vpow(:,1).*w1.^vpow(:,2));
% Integral of polynomial P(x1,wq) with respect to probability measure of w1
Pcc=Co'*(x1.^vpow(:,1).*yw1(vpow(:,2)+1));

%% Plots

% plot W
[x1,w1]=meshgrid([-0.9:0.01:0.9],[-0.9:0.01:0.9]);
surf(x1,w1,eval(P),'FaceColor','red','FaceAlpha',0.85,'EdgeColor','none','FaceLighting','phong');hold on;grid on;
surf(x1,w1,ones(size(eval(Pcc))),'FaceColor','blue','FaceAlpha',0.5,'EdgeColor','none','FaceLighting','phong')
axis tight;hidden off;ylabel('w1');xlabel('x1');hold on;xlim([-1,1]);ylim([-1,1]);zlim([0,3]);
camlight; lighting gouraud,hold on
title('$P(x_1,\omega_1)$','Interpreter','latex', 'FontSize',31);
zlim([0 2])

% plot set K1
for x=0:0.05:1; for q=0:0.05:1
        if (0.5*q*(q^2+(x-0.5)^2)-(q^4+q^2*(x-0.5)^2+(x-0.5)^4)) >=0; plot3([x x],[q q],[0 1],'k-*');end
end;end
xlabel('$x_1$','Interpreter','latex', 'FontSize',31);ylabel('$\omega_1$','Interpreter','latex', 'FontSize',31)
str1 = '$ \mathcal{K} $';text(0.94,0.2,str1,'HorizontalAlignment','right','Interpreter','latex','FontSize',30) 

% plot Integral of W 
figure
x1=[-1:0.01:1];
plot(x1,eval(Pcc),x1,1-Delta*ones(size(x1)),'LineWidth',5);grid;hold on
%title('$ \int \mathcal{P}^d_{\mathcal{W}}(x,q) d\mu_q$, \ \ \beta$','Interpreter','latex', 'FontSize',31);
xlabel('$x_1$','Interpreter','latex', 'FontSize',31);set(gca,'fontsize',20)
str2 = '$ 1-\Delta $';text(0.5,0.5,str2,'HorizontalAlignment','right','Interpreter','latex','FontSize',30)         
str3 = '$ \int {\mathcal{P}}(x_1,\omega_1) d\mu_{\omega_1}$';text(-0.3,0.4,str3,'HorizontalAlignment','right','Interpreter','latex','FontSize',30)  

ylim([0 1])

