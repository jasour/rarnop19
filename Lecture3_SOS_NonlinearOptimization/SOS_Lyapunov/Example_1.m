% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 3: Sum Of Squares Based SDP For Nonlinear Optimization
%% SOS Lyapunov Function Search for nonliner system dx(t)/dt=f(x):   V(0)=0, V(x):SOS  -dV(x)/dt:SOS
clc;clear all;close all

% variables
sdpvar x1 x2
% nonlinear system dx(t)/dt=f(x)
f=[-x1+(1+x1)*x2;-(1+x1)*x1];
x = [x1;x2];
% Polynomial Lyapunov Function or order 4 
[V,c] = polynomial(x,4);
% dV(x)/dt
dVdt = jacobian(V,x)*f;
% SOS COnditions
F = [sos(V),sos(-dVdt),c(1)==0];
%SDP solver
ops = sdpsettings('solver','mosek');
% SOlve SOS program
[sol,v,Q]=solvesos(F,[],ops,c);


%% Results

% obtained PSD matrix Q1: v(x)=B(x)'Q1 B(x)
Q1=value(Q{1});
% obtained basis B(x):  v(x)=B(x)'Q1 B(x)
B1=sdisplay(v{1});
% obtained PSD matrix Q2: -dv(x)/dt=B(x)'Q2 B(x)
Q2=value(Q{2});
% obtained basis B(x):  -dv(x)/dt=B(x)'Q2 B(x)
B2=sdisplay(v{2});

% obtains sum of squares polynomials: v(x)=h(x)'h(x)
h=sosd(F); sdisplay(h);
hh=sdisplay(sosd(F(1))'*sosd(F(1)));

% plot Lyapunov Function v(x) and dv(x)/dt
L1=sdisplay(v{1}'*Q{1}*v{1});
dL1=sdisplay(v{2}'*Q{2}*v{2});

L2=strrep(strrep(L1,'*','.*'),'^','.^');L3=cell2mat((L2));
dL2=strrep(strrep(dL1,'*','.*'),'^','.^');dL3=cell2mat((dL2));

[x1,x2]=meshgrid([-1:0.01:1],[-1:0.01:1]);
surf(x1,x2,eval(L3),'FaceColor','red','FaceAlpha',0.85,'EdgeColor','none','FaceLighting','phong');hold on;grid on;
surf(x1,x2,-eval(dL3),'FaceColor','blue','FaceAlpha',0.85,'EdgeColor','none','FaceLighting','phong');hold on;grid on;

camlight; lighting gouraud