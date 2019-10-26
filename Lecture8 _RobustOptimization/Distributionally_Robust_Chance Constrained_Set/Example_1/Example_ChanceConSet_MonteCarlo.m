
Prob_List=[];
for x=-1:0.01:1

% number of uncertainty samples
N=1000000;

q=random('Normal',0,0.2,1,N);

K=0.8^2-x.^2-q.^2;
Pro=size(find(K>=0),2)/N;
Prob_List=[Prob_List;x,Pro];

end
plot([-1:0.01:1],Prob_List(:,2),'color','blue','LineWidth',3)
plot([-1:0.01:1],1-Delta*ones(size(x1)),'--','color','blue','LineWidth',3)
str2 = '$ 1-\Delta $';text(0.5,0.8,str2,'HorizontalAlignment','right','Interpreter','latex','FontSize',30,'color','blue')         
hold on; xlim([-0.9,0.9]);ylim([0 1.5])