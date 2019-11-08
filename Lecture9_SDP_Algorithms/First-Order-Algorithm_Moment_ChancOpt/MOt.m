function out = MOt(y)

global AInd_y
n = sqrt(length(y));
X = reshape(y, n, n);
out=AdjMoment(X,AInd_y);