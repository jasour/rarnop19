
Prob_List=[];
for x=-1:0.01:1

% number of uncertainty samples
N=1000000;

q=random('Uniform',-1,1,1,N);

K=(-0.0001+0.0032*x+0.0027*q-0.0384*x^2-0.0440*x.*q-0.0030*q.^2+0.2048*x.^3+0.1760*x.^2*q+0.0484*x.*q.^2+0.0832*q.^3-0.4096*x.^4-0.1936*x.^2.*q.^2-0.0915*q.^4);

Pro=size(find(K>=0),2)/N;
Prob_List=[Prob_List;x,Pro];

end
plot([-1:0.01:1],Prob_List(:,2),'--','LineWidth',3)

 ylim([0 1])