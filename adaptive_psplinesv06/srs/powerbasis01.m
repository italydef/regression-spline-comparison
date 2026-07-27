function xm = powerbasis01(x,degree,knots,der) ;
%
%	Returns the power basis functions of a spline of given degree
%	USAGE: xm = powerbasis(x,degree,knots) 
%
%	Last edit: 	3/22/99
%

if nargin < 4 ;
	der = 0 ;
end ;

if der > degree ;
disp('********************************************************') ;
disp('********************************************************') ;
disp('WARNING:  der > degree --- xm not returned by powerbasis') ;
disp('********************************************************') ;
disp('********************************************************') ;
return ;
end ;

n=size(x,1) ;
nknots = length(knots);
mx = mean(x) ;
stdx = std(x) + 100*eps ;

if der == 0 ;
	xm=ones(n,1);
	else ;
	xm = zeros(n,1) ;
end ;

	for i=1:degree ;
		if i < der ;
			xm = [xm zeros(n,1)] ;
			else ;
			xm = [xm prod((i-der+1):i) *  x.^(i-der)] ;
		end ;
	end ;
	
	if nknots > 0 ;
		for i=1:(nknots) ;
		xm = [xm prod((degree-der+1):degree) * ...
			 (x-knots(i)).^(degree-der).*(x > knots(i))] ;
		end ;
	end ;
