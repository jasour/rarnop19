% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 8: Nonlinear Robust Optimization
%% Monte-Carlo Test of Example 2

clc;clear;


X_R_m=[];
for u=-1:0.01:1

% number of uncertainty samples
N=100000;

% s is (n x N), N points in n-ball with radius R
n = 4;R=0.1;
s = randn(n,N);
r = rand(1,N).^(1/n);
c = r./sqrt(sum(s.^2,1));
s = bsxfun(@times, s, c)*R;

x1=s(1,:);x2=s(2,:);x3=s(3,:);w=s(4,:);

% x(k+1)
f1=0.2*w.*x2;
f2=x1.*x3;
f3=1.2*x1-0.5*x2+x3+2*u;

% control objective
g=-(1^2-(f1-0).^2/0.2^2-(f2-0).^2/0.2^2-(f3-0.9).^2/0.4^2);

Pro=size(find(g<=0),2)/N;
if Pro==1; X_R_m=[X_R_m,u]; end

end



plot(X_R_m,zeros(size(X_R_m)),'r-','LineWidth',3);
xlim([-2 2]); xlabel('$u$','Interpreter','latex', 'FontSize',40);
str1 = '$U_{R}$';text(0.1,0.1,str1,'Interpreter','latex','FontSize',30,'color','red')
grid on;xlim([0 1]);ylim([-0.5 1]);


