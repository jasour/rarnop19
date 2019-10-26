function y = glex2num(x)
% function y = glex2num(x)
% X = [X_1,...,X_n] are non-negative integers
% Y returns the index number of X in graded lexicographic order

% Chao Feng, May 10, 2009
% Last modified: July 15, 2010

n = length(x);
glex = sum(x);

if glex == 0,
    y = 1;
else
    y = nchoosek(glex-1+n,n) + glex2num(x(2:end));
end
