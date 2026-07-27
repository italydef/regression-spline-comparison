function y = cdfpois(x,mu) ;
%  Compute the cdf at x of a Poisson random variable with mean mu.
%  x 
n = rows(x) ;
m = cols(x) ;
y = zeros(n,m) ;
x = floor(x) ;
	for i = 1:n ;
		for j = 1:m
		if x(i,j) < 0 ;
		  y(i,j) = 0 ;
		else
		  y(i,j) = sum(pdfpois([0:x(i,j)],mu)) ;
		end ;
		end ;
	end ;




