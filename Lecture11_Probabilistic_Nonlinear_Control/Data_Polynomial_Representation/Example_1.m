% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 11: Risk Aware Planning and Control Of Probabilistic Nonlinear Dynamical Systems
%% Semialgebraic Set Representation of Data

clc;clear;close all; 


% order of polynomial 
d=20;
d_sos=d;



% variables
sdpvar x1 x2

% polynomial  with unknown coefficients coef
[P,co_p] = polynomial([x1 x2],d);

% Simple Box
B1=(x1+1)*(1-x1); 
B2=(x2+1)*(1-x2);


% Random Data
s=100;
xs1=random('Normal',0.5,.1,1,s);
ys1=random('Normal',0.3,.1,1,s);

xs2=random('Normal',-0.3,.1,1,s);
ys2=random('Normal',-0.5,.1,1,s);

xs3=random('Normal',-0.3,.1,1,s);
ys3=random('Normal',0.4,.1,1,s);

xs=[xs1,xs2,xs3];ys=[ys1,ys2,ys3];
plot(xs,ys,'.');hold on

% Linear Constraints
F=[];pow=[];for i=0:d; pow=[pow;genpow(2,i)];end
for i=1:size(xs,2)
F=[F, sum(co_p.*(xs(i).^pow(:,1)).*(ys(i).^pow(:,2)))>=1];
end

% SOS Constraints
[s1,c1] = polynomial([x1 x2],d_sos);
[s2,c2] = polynomial([x1 x2],d_sos);
F = [F, sos(P-[s1]*B1-[s2]*B2), sos(s1), sos(s2)];

% Integral
Int=int(P,[x1 x2],[-1 -1],[1 1]);


ops = sdpsettings('solver','mosek');
[sol,v,Q]=solvesos(F, Int,[],[c1;c2;co_p]);


%% Results
% coefficients of polynomial
C=value(co_p');

% Polynomial
syms x1 x2
PP= sum(C'.*(x1).^pow(:,1).*(x2).^pow(:,2));
[x1,x2]=meshgrid([-0.99:0.01:0.99],[-0.99:0.01:0.99]);
surf(x1,x2,eval(PP),'FaceColor','blue','FaceAlpha',0.5,'EdgeColor','none','FaceLighting','phong');hold on
camlight; lighting gouraud

[x1,x2]=meshgrid([-0.99:0.01:0.99],[-0.99:0.01:0.99]);
surf(x1,x2,ones(size(eval(PP))),'FaceColor','red','FaceAlpha',0.5,'EdgeColor','none','FaceLighting','phong');hold on
contour(x1,x2,eval(PP),'ShowText','on')

figure(2);plot(xs,ys,'.');hold on
contour(x1,x2,eval(PP),[1 1],'--rs',...
    'LineWidth',2)
