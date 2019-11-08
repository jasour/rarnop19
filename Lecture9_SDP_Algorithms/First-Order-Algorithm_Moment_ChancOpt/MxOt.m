function out = MxOt(y)

global AInd_x
n = sqrt(length(y));
X = reshape(y, n, n);
out= AdjMoment(X,AInd_x);