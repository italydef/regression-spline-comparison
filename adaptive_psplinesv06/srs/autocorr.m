function y = autocorr(x,k) ;
% Computes the first k autocorrelation coefficients of x
% 	Usage: y = autocorr(x,k) ;
%
%		INPUT (required)
%	x = time series
%
%		INPUT (optional)
%	k = number of lags at which the autocorr is calculated
%		(default = 20)
if nargin < 2 ;
	k = 20 ;
end ;

if size(x,1) == 1 ;
	x = x' ;
end ;

x2 = x ;
y = zeros(k,1) ;
  for i = 1:k ;
  n = rows(x) ;
  x = x(1:(n-1)) ;
  x2 = x2(2:n) ;
  c = corrcoef(x,x2) ;
  y(i) = c(1,2) ;
  end ;
