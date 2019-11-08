function out = UPOt(y)

global Yq Mind_xq AInd_y Nx
n = sqrt(length(y));
X = reshape(y, n, n);
out = AdjMoment_xq(X,Yq,Mind_xq,AInd_y,Nx);