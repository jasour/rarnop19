% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 8: Nonlinear Robust Optimization

clc;clear;

X_R_m=[];
for x=-1:0.01:1
N=1000000;
q=random('Uniform',-0.5,0.5,1,N);
K=x.^2+q.^2-1;
Pro=size(find(K<=0),2)/N;
if Pro==1; X_R_m=[X_R_m,x]; end
end


figure;
plot(X_R_m,zeros(size(X_R_m)),'--','LineWidth',3);
xlim([-2 2]); xlabel('$x$','Interpreter','latex', 'FontSize',40);
str1 = '$\chi^d_R$';text(0,0.3,str1,'Interpreter','latex','FontSize',30,'color','black')
grid on


