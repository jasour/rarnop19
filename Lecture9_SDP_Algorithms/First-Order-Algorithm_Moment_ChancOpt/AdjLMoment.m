function AdjLM=AdjLMoment(X,AInd_L)
global Ny Lmo Lco
% AdjLM=zeros(size(AInd_L,1),1);
% for i=1:size(AInd,1)
%     AdjM(i)= sum(X(AInd(i,1:nnz(AInd(i,:)))));
% end
% 

AdjLM=zeros(Ny,1);
for i=1:size(Lmo,3); YY(:,:,i)=X; end
YY=YY.*Lco;
for i=1:Ny
%[ind]=find(Lmo==i);
%if isempty(ind)~=1; Adj(i)=sum(YY(ind)); end
AdjLM(i)=sum(YY(AInd_L(i,1:nnz(AInd_L(i,:)))));
end
