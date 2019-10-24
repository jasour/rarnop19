% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 8: Nonlinear Robust Optimization
%% SOS Barrier Function based Safety Verification of "uncertain" nonliner system dx(t)/dt=f(x,u,w) with respect to unsafe set X_obs 
% w is a bounded uncertainty, u is a given control policy. 
%%
clc;clear all;close all

% UAV Dynamics
% states=[x,y,phi]: [position, yaw angle];  
% Dynamics: f1=v*sin(phi)+omega; f2=v*cos(phi); f3=u;
% v=1 m/s is the speed of the airplane
% omega: uncertainty: “cross-wind” (bounded between [-0.05, 0.05]). 
% control policy: states feedback : u=-K*(phi-phi_des) where phi_des is
% desired direction 

% tf: execution time of control policy
% y_0: initial states [x,y,phi]
tf=5; y0 = [0;0;0]; v=1; phi_des= 0*pi/180;
tspan = [0 tf]; [t,y] = ode45(@(t,y) [v*sin(y(3))+0.2*rand-0.1;v*cos(y(3));-50*(y(3)-phi_des)], tspan, y0);
plot(y(:,1),y(:,2),'LineWidth',2);xlabel('x');ylabel('y');
hold on; ax = gca; ax.FontSize = 30;
str1 = '$(x_t,y_t,\phi_t)$';text(0,0,str1,'Interpreter','latex','FontSize',30)


% Polynomial dynamics: Taylor expansion of nonlinear dynamics to degree 3  
syms x y phi omega
f1=v*sin(phi)+omega;
f2=v*cos(phi);
u=-50*(phi-phi_des); f3=u;
% Obtained Polynomial dynamics
fT1 = taylor(f1, phi, 'Order', 3);
fT2 = taylor(f2, phi, 'Order', 3);
fT3=f3;


%% SOS Program

% 2d: Order of polynomial Barrier function
d=2;

% variables
sdpvar x y phi omega

% Polynomial Uncertain nonlinear system dx(t)/dt=f(x,w)
f=[eval(fT1);eval(fT2);eval(fT3)];

% polynomial Barrier function of order 2d, c: coefficients, Vm: monomials
[V,c,Vm] = polynomial([x;y;phi],2*d);

% dV(x)/dt
dVdt = jacobian(V,[x;y;phi])*f;

% Uncertainty Set [-0.1 0.1]
g_omega=(0.1-omega)*(omega+0.1);
% Obstacle set 
g_obs=0.25^2-(x-0.5)^2-(y-0.5)^2;
% Set of all states 
X=2^2-(x)^2-(y)^2-(phi)^2;



% s1: SOS polynomial
[s1,c1] = polynomial([x;y;phi;omega],2*d);
% s2: SOS polynomial
[s2,c2] = polynomial([x;y;phi;omega],2*d);
% s3: SOS polynomial
[s3,c3] = polynomial([x;y;phi;omega],2*d);



% SOS Conditions:
%V([x_0,y_0,phi_0])=0
% V>=1 on Obstacle  set 
% -dVdt>=0 for all uncertainty w  
F = [sos(V-1-s1*g_obs),sos(-dVdt-s2*g_omega-s3*X), sos(s1),sos(s2),sos(s3),c(1)==0];

%SDP solver
ops = sdpsettings('solver','mosek');

% Solve SOS based SDP
[sol,v,Q]=solvesos(F,[],ops,[c1;c2;c3;c]);
%%

%% Results
% Obtained coefficients of polynomial Lyapunov Function
cc=value(c);

% Lyapunov Function
L1=sdisplay(cc'*Vm);
% dV(x)/dt
dL1 = sdisplay(jacobian(cc'*Vm,[x;y;phi])*f);

%% Plots

L2=strrep(strrep(L1,'*','.*'),'^','.^');L3=cell2mat((L2));
dL2=strrep(strrep(dL1,'*','.*'),'^','.^');dL3=cell2mat((dL2));

[x,y,phi]=meshgrid([-1:0.01:1],[-1:0.01:1],[-1:0.01:1]);

p=patch(isosurface(x,y,phi,eval(L3),0));
 set(p,'FaceColor','red','EdgeColor','none','FaceAlpha',0.7);

p=patch(isosurface(x,y,phi,eval(L3),1));
set(p,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.7);
camlight; lighting gouraud

ezplot('0.25^2-(x-0.5)^2-(y-0.5)^2');

plot(0,0,'ko','LineWidth',2)

xlim([-1 1]);ylim([-1 1]);title('Barrier Function')


xlabel('$x$','Interpreter','latex', 'FontSize',40);ylabel('$y$','Interpreter','latex', 'FontSize',40);zlabel('$\phi$','Interpreter','latex', 'FontSize',40);
str1 = '$V(x)\leq 0$';text(-0.5,0.5,0.9,str1,'Interpreter','latex','FontSize',30)
str1 = '$V(x)\geq 1$';text(0.5,-0.5,0.5,str1,'Interpreter','latex','FontSize',30)
str1 = '$\chi_{obs}$';text(0.5,0.5,0,str1,'Interpreter','latex','FontSize',30)
 
[x,y]=meshgrid([-1:0.01:1],[-1:0.01:1]);
surf(x,y,zeros(size(x)),'FaceColor','black','FaceAlpha',0.4,'EdgeColor','none','FaceLighting','phong');hold on;grid on;


