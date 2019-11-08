function AdjM=AdjMoment(X,AInd)

AdjM=zeros(size(AInd,1),1);
for i=1:size(AInd,1)
    AdjM(i)= sum(X(AInd(i,1:nnz(AInd(i,:)))));
end
