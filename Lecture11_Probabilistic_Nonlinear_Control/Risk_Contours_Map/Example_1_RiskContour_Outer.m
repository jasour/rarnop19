% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 11: Risk Aware Planning and Control Of Probabilistic Nonlinear Dynamical Systems
%% Risk Contours Map
% To obtain the risk contours map described in page 52 of Lecture 11, we
% solve dual Optimization of SOS-SDP in page 160, Lecture 7:Nonlinear Chance Constrained and Chance Optimization.
% Dual read as page 107, Lecture 10: Probabilistic Nonlinear Safety Verification.
clc;clear all; close all
mset clear; mset('yalmip',true); mset(sdpsettings('solver','mosek'))

% polynomial order
d=10;
% polynomial power
vpow=[];for k = 0:2*d; vpow = [vpow;genpow(3,k)]; end

% risk levels
Delta = [0.9, 0.91, 0.93, 0.95 0.97 0.99];



%%
% location x and y, uncertainty w
mpol x y w
% assigned measure
mu = meas(x,y,w);
% moments
m = mom(x.^vpow(:,1).*y.^vpow(:,2).*w.^vpow(:,3));

%slack variables
mpol xs ys ws
% assigned measure
mus = meas(xs,ys,ws); 
% moments
ms = mom(xs.^vpow(:,1).*ys.^vpow(:,2).*ws.^vpow(:,3));

% moments of Lebesgue measure over [-1,1]     
m_leb=[2];for i=1:2*d ;m_leb(i+1,1)=(1/1)*(((1)^(i+1) - (-1)^(i+1))/(i+1));end

% moment of uncertainty w \in [0,1]: Beta distribution with parameters alpha and beta
alpha=1.1;beta=5;
m_w=[1];for k=1:2*d; m_w=[m_w,(alpha+k-1)/(alpha+beta+k-1)*m_w(end) ]; end
m_w=m_w';
xxx = 0:.001:1;yyy = betapdf(xxx,alpha,beta);plot(xxx,yyy);title('Beta distribution')

%upper bound moments: moments of (Lebesgue measure on x: [-1 1])*(Lebesgue measure
%on y [-1 1])*(uncertainty w)
m_up=m_leb(vpow(:,1)+1).*m_leb(vpow(:,2)+1).*m_w(vpow(:,3)+1);
     
% Uncertain Safe set
a=2.5;a2=1.5;sx=1;sy=1; 
 Xsafe=[-(-1*(a*(sx*x))^4+0.5*((a*(sx*x))^2-(a2*((sy*y)+0)).^2)+0.01+0.5*w); (3-x^2-y^2-w^2) ]      
   
 % moment SDP
 P = msdp(max(mass(mu)),m+ms == m_up, Xsafe>=0, w*(1-w)>=0, (1-x^2)>=0, (1-y^2)>=0,(1-xs^2)>=0, (1-ys^2)>=0, ws*(1-ws)>=0);
 
 % solve moment SDP
[stat,obj,mm,dual] = msol(P);

%% Reults

if stat >=0
% obtained dual variables: coefficients of polynomial in x,y,w
dual;

% polynomial construction
syms x y
PP=sum(dual.*x.^vpow(:,1).*y.^vpow(:,2).*m_w(vpow(:,3)+1));
[x,y] = meshgrid(-0.9:0.01:0.9,-0.9:0.01:0.9);
PP=eval(PP);
surf(x,y,PP); axis square
xlim([-0.9 0.9]);ylim([-0.9 0.9]);zlim([0.7 1.1])

end
%% Plots
figure; hold on

        
        F=contour(x,y,PP,Delta(1)*[1,1],'r','LineWidth',1,'ShowText','off'); hold on;
        if size(F,2)>0 ;F=[F(:,ceil(end/2):end)]; plot3(F(1,:),F(2,:),Delta(1)*ones(1,size(F,2)),'y','LineWidth',2); end
        
        F=contour(x,y,PP,Delta(2)*[1,1],'r','LineWidth',1,'ShowText','off'); hold on;
        if size(F,2)>0 ;F=[F(:,ceil(end/2):end)]; plot3(F(1,:),F(2,:),Delta(2)*ones(1,size(F,2)),'y','LineWidth',2); end

        
        F=contour(x,y,PP,Delta(3)*[1,1],'r','LineWidth',1,'ShowText','off'); hold on;
        if size(F,2)>0 ;F=[F(:,ceil(end/2):end)]; plot3(F(1,:),F(2,:),Delta(3)*ones(1,size(F,2)),'y','LineWidth',2); end
   
        F=contour(x,y,PP,Delta(4)*[1,1],'r','LineWidth',1,'ShowText','off'); hold on;
        if size(F,2)>0 ;F=[F(:,ceil(end/2):end)]; plot3(F(1,:),F(2,:),Delta(4)*ones(1,size(F,2)),'y','LineWidth',2); end
        
surf(x,y,PP,'FaceColor','red','EdgeColor','none','FaceAlpha',0.7);

camlight; lighting gouraud
xlim([-0.9 0.9]);ylim([-0.9 0.9]);zlim([0.7 1.1])
hold on;grid on;set(gca,'fontsize',31)
xlabel('$x_1$','Interpreter','latex', 'FontSize',31)
ylabel('$x_2$','Interpreter','latex', 'FontSize',31)
zlabel('$\mathcal{P}_{inner}(x)$','Interpreter','latex', 'FontSize',31)
axis square


figure
        contour(x,y,PP,Delta(1)*[1,1],'k','LineWidth',1,'ShowText','off'); hold on;
        contour(x,y,PP,Delta(2)*[1,1],'--k','LineWidth',1,'ShowText','off'); hold on;
        contour(x,y,PP,Delta(3)*[1,1],'k','LineWidth',1,'ShowText','off'); hold on;
        contour(x,y,PP,Delta(4)*[1,1],'--k','LineWidth',1,'ShowText','off'); hold on;
xlim([-0.9 0.9]);ylim([-0.9 0.9]);
set(gca,'fontsize',31)
xlabel('$x_1$','Interpreter','latex', 'FontSize',31)
ylabel('$x_2$','Interpreter','latex', 'FontSize',31)
axis square

figure;
Fs=20; Fs1=12; Fst=20;
l=-0.9;u=0.9;
subplot(2,3,1);hold on    
        contour(x,y,PP,Delta(1)*[1,1],'r--','LineWidth',2);xlim([l u]);ylim([l u]);set(gca,'fontsize',Fs1)
        xlabel('$x_1$','Interpreter','latex', 'FontSize',Fs);ylabel('$x_2$','Interpreter','latex', 'FontSize',Fs)
        title('$\bar{\mathcal{C}}^{\Delta=0.05}_{r}, \ {\mathcal{C}}^{\Delta=0.05}_{r}$','Interpreter','latex', 'FontSize',Fst)
        axis square
subplot(2,3,2);hold on
contour(x,y,PP,Delta(2)*[1,1],'r--','LineWidth',2);xlim([l u]);ylim([l u]);set(gca,'fontsize',Fs1)
        xlabel('$x_1$','Interpreter','latex', 'FontSize',Fs);ylabel('$x_2$','Interpreter','latex', 'FontSize',Fs)
        title('$\bar{\mathcal{C}}^{\Delta=0.1}_{r}, \ {\mathcal{C}}^{\Delta=0.1}_{r}$','Interpreter','latex', 'FontSize',Fst)
        axis square

subplot(2,3,3);hold on
        contour(x,y,PP,Delta(3)*[1,1],'r--','LineWidth',2);xlim([l u]);ylim([l u]);set(gca,'fontsize',Fs1)
        xlabel('$x_1$','Interpreter','latex', 'FontSize',Fs);ylabel('$x_2$','Interpreter','latex', 'FontSize',Fs)
        title('$\bar{\mathcal{C}}^{\Delta=0.2}_{r}, \ {\mathcal{C}}^{\Delta=0.2}_{r}$','Interpreter','latex', 'FontSize',Fst)
        axis square
subplot(2,3,4);hold on
        contour(x,y,PP,Delta(4)*[1,1],'r--','LineWidth',2);xlim([l u]);ylim([l u]);set(gca,'fontsize',Fs1)
        xlabel('$x_1$','Interpreter','latex', 'FontSize',Fs);ylabel('$x_2$','Interpreter','latex', 'FontSize',Fs)
        title('$\bar{\mathcal{C}}^{\Delta=0.3}_{r}, \ {\mathcal{C}}^{\Delta=0.3}_{r}$','Interpreter','latex', 'FontSize',Fst)
        axis square
