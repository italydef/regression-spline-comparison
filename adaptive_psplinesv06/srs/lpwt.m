function y = lpwt(x,p) ;
%
%  weights for l-p norm penalty
%
y = (abs(x) + 100*eps).^(p-2) ;
