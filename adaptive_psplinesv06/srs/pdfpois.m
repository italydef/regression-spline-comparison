function y = pdfpois(x,mu) ;
%  Compute the pdf at x of a Poisson random variable with mean mu.
%
if mu == 0 ;
  y = (x == 0) ;
else ;
  n = rows(x) ;
  m = cols(x) ;
  y = zeros(n,m) ;
	for i = 1:n ;
		for j = 1:m
		if x(i,j) < 0 ;
		  'WARNING: x is negative' 
		  y(i,j) = 0 ;
		  elseif x(i,j) ~= ceil(x(i,j));
		  'WARNING: x is not an integer'
		  y(i,j) = 0;
		  else;
		  x2=max([x(i,j) 1]);
		  y(i,j) =  exp( -mu + log(mu)*x(i,j) -   ...
		  sum(log([1:x2])) ) ;
		  end ;
		end ;
	end ;
end ;



