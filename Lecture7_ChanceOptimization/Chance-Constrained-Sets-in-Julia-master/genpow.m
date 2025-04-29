function v = genpow(ndig,sum)
% GENPOW - Internal use only
  
% Generate powers
%
% GENPOW(NDIG,SUM) returns a matrix whose rows are
% all vectors with NDIG ndigits summing up to SUM
% For example GENPOW(3,2) returns [2 0 0;1 1 0;1 0 1;0 2 0;0 1 1;0 0 2]
%
% Used by GENIND

% D. Henrion, 20 November 2003
% Last modified on 16 December 2006
  
if ndig > 1
 v = zeros(1,ndig);
 if sum > 0
  r = 0;
  for k = sum:-1:0
   % Recursive call
   w = genpow(ndig-1,sum-k); rd = size(w,1);
   v(r+1:r+rd,1) = repmat(k,rd,1);
   v(r+1:r+rd,2:ndig) = w;
   r = r+rd;
  end
 end
else
 v = sum;   
end

