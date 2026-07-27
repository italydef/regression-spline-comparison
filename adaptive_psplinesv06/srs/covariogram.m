function [d,c,v]=covariogram(x,y) ;
%
%  USAGE:  [d,c,v]=covariogram(x,y) ;
%   d = distances
%   c = covariances
%   v = variogram
%
n=rows(x) ;
mu = mean(y) ;
d = zeros(n*(n-1)/2,1) ;
c = d ;
v = c ;
	k= 1 ;
	for i = 1:(n-1) ;
		for j=(i+1):n ;
		d(k) = sqrt( sum((x(i,:)-x(j,:)).^2) );
		c(k) = (y(i)-mu)*(y(j)-mu) ; 
		v(k) = (y(i) - y(j)).^2 ;
		k=k+1;
		end ;
	end ;
