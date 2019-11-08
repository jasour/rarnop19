function out = LOt(y)

global AInd_L
n = sqrt(length(y));
LMc = reshape(y, n, n);
%out= AdjLM(LMc,LMInd1);

out=AdjLMoment(LMc,AInd_L); 