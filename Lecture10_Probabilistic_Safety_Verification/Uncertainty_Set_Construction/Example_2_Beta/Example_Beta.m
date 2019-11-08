%% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 10: Probabilistic Safety Verification
%% Uncertainty Set Construction from the moment information of the uncertainty
% Uncertainty Set: K_d={x:  p_d(x) > = 1}, d-th order polynomial p_d(x) is obtained by
% solving a SOS program
% Ashkan Jasour, C. Lagoa, ”Reconstruction of Support of a Measure From Its Moments”, 53 rd IEEE Conference on Decision and Control, Los Angeles, California, 2014
%%
clc;clear;close all;


% Order of polynomial
n=50;

%% Parameters
% dimension of x
nx= 1;
% Order of Localization Matrix
d=1*ceil(n)+0;
% support of uniform measure
b2=1;b1=0;
% Simple Box
B2=1.2;B1=-1.2;
w=0;
%% Varibales
x=sdpvar(nx);
sdpvar I
v0 = monolist([x],n/2); Q0 = sdpvar(length(v0));
v1 = monolist([x],(n-2)/2);Q1 = sdpvar(length(v1));
p_sos = v0'*Q0*v0 + (x-B1)*(B2-x)*v1'*Q1*v1; % possitive on the Simple Box

[f,a] = polynomial([x],n); % indictor 


%% Moments of Beta dist
Alpha = 4;Beta= 4;
yq=[1];for i=1:2*d ;yq(i+1,1)= b2*((Alpha+i-1)/(Alpha+Beta+i-1))*yq(end);end ; m= yq;

%% Cost Function
J=int(f,x,B1,B2);

%% Localization Matrix
d1=floor(abs(2*d- n )/2);LM=0;
for i=1:n+1
H = hankel([i:i+(d1)],[i+d1:i+2*(d1)]);
if i==1; LM=LM+(a(i)-I)*m(H); else LM=LM+a(i)*m(H); end
end
%% SDP
F = [coefficients(f-p_sos,[x]) == 0, Q0 >= 0, Q1 >= 0, LM >=0, I>=1, I<1.2]; 
solvesdp(F,J-1.2*I,[]);

%% Solution
double(Q0);
double(J)
aa=double(a);

x=[-3:0.01:3];
I=0;
for i=1:n+1
I=I+aa(i)*x.^(i-1);
end
plot(x,I,'m','LineWidth',2);xlim([B1,B2]);ylim([-0.1,2.5]);hold on;grid on; %axis square
xlabel('x');% ylabel('')
plot([B1,b1,b1,b2,b2,B2],[0 0 1 1 0 0],'--k','LineWidth',2);%plot([b1,b2],[1 1],'ro')

