function out = LO(x)

global Lmo Lco
%LM = LMoment(Dc1,Cc1,x,d1);
LM=sum(Lco(:,:,:).*x(Lmo(:,:,:)),3);
out = LM(:);
