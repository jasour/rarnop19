% MIT 16.S498: Risk Aware and Robust Nonlinear Planning, Fall 2019
% Lecture 8: Nonlinear Robust Optimization
%%
clc;

Prob_List=[];
for x=-1:0.01:1

% number of uncertainty samples
N=100000;

% random mean \in [-0.1,0.1] and random sigma \in [0.1,0.3] 
m=random('Uniform',-0.1,0.1,1,N);
sigma=random('Uniform',0.1,0.3,1,N);
% random q with Normal probability distribution with uncertain mean and
% sigma
q=random('Normal',m,sigma);

% Safety Constraint K={g(x,q)>=0}
K=0.8^2-x.^2-q.^2;

% Probability of design variable x 
Pro=size(find(K>=0),2)/N;
Prob_List=[Prob_List;x,Pro];

end

plot([-1:0.01:1],Prob_List(:,2),'color','blue','LineWidth',3);hold on
% Distrbutionally Robust Chance Constrained Set 
% X_DR = {all x:   probability { g(x,q)>=0 }>=1-Delta, for all probability distribution of uncertainty parameter q } 
% Probability >= 1-Delta
plot([-1:0.01:1],1-Delta*ones(size(x1)),'--','color','blue','LineWidth',3)
str2 = '$ 1-\Delta $';text(0.4,0.95,str2,'Interpreter','latex','FontSize',30,'color','blue')         
str3 = '$ Prob Curve (Monte-Carlo)$';text(-0.2,1.1,str3,'Interpreter','latex','FontSize',30,'color','blue')  
ylim([0 1.5]); grid on