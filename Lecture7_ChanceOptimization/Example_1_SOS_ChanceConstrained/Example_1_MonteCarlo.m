
Prob_List=[];
for x=-1:0.01:1

% number of uncertainty samples
N=1000000;

q=random('Uniform',-1,1,1,N);

K=0.5*q.*(q.^2+(x-0.5).^2)-(q.^4+q.^2.*(x-0.5).^2+(x-0.5).^4);
Pro=size(find(K>=0),2)/N;
Prob_List=[Prob_List;x,Pro];

end
plot([-1:0.01:1],Prob_List(:,2),'--','LineWidth',3)

 ylim([0 1])