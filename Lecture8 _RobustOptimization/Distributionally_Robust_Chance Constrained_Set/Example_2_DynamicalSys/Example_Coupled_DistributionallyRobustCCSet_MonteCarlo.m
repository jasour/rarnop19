% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 8: Nonlinear Robust Optimizaion
%% Distributionally Robust Chance Constrained Set 
%% Monte Carlo test

Prob_List=[];
for u=-1:0.01:1

% number of uncertainty samples
N=100000;

n = 4;R=0.1;
s = randn(n,N);
r = rand(1,N).^(1/n);
c = r./sqrt(sum(s.^2,1));
s = bsxfun(@times, s, c)*R;
m1=s(1,:);m2=s(2,:);m3=s(3,:);m4=s(4,:);

n = 4;R=0.1;
s = randn(n,N);
r = rand(1,N).^(1/n);
c = r./sqrt(sum(s.^2,1));
s = bsxfun(@times, s, c)*R;
sigma1=s(1,:)+0.2;sigma2=s(2,:)+0.2;sigma3=s(3,:)+0.2;sigma4=s(4,:)+0.2;


% m1=random('Uniform',-0.1,0.1,1,N);sigma1=random('Uniform',0.1,0.3,1,N);
% m2=random('Uniform',-0.1,0.1,1,N);sigma2=random('Uniform',0.1,0.3,1,N);
% m3=random('Uniform',-0.1,0.1,1,N);sigma3=random('Uniform',0.1,0.3,1,N);
% m4=random('Uniform',-0.1,0.1,1,N);sigma4=random('Uniform',0.1,0.3,1,N);

x1=random('Normal',m1,sigma1);
x2=random('Normal',m2,sigma2);
x3=random('Normal',m3,sigma3);
w=random('Normal',m4,sigma4);

% x(k+1)
f1=0.2*w.*x2;
f2=x1.*x3;
f3=1.2*x1-0.5*x2+x3+2*u;

% success set
K=-(1^2-(f1-0.1).^2/0.2^2-(f2-0.1).^2/0.2^2-(f3-0.4).^2/0.3^2);


Pro=size(find(K>=0),2)/N;
Prob_List=[Prob_List;u,Pro];

end
plot([-1:0.01:1],Prob_List(:,2),'--','LineWidth',3,'color','blue');hold on
plot([-1:0.01:1],1-Delta*ones(size([-1:0.01:1])),'--','color','blue','LineWidth',5);grid;hold on
str2 = '$ 1-\Delta $';text(0.6,0.9,str2,'Interpreter','latex','FontSize',30,'color','blue')         
grid on; ylim([0 1.5])