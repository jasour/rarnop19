
Prob_List=[];
for x=-1:0.01:1

% number of uncertainty samples
N=1000000;

q=random('Uniform',-1,1,1,N);

K=-[((3*x)/2 + 1).^2/4 - ((3*x)/2 + 1).^3/4 + ((3*x)/2 + 1).^4/16 + (9*q.^2)/100 - 29/400];
Pro=size(find(K>=0),2)/N;
Prob_List=[Prob_List;x,Pro];

end
plot([-1:0.01:1],Prob_List(:,2),'--','LineWidth',3)

 ylim([0 1])