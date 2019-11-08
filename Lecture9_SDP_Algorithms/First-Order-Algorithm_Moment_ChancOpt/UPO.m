function out = UPO(x)

global Yx Yq
MAxT = MomentCross(x,Yx,Yq); 
out = MAxT(:);



