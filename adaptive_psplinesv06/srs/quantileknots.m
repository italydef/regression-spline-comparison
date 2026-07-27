function knots = quantileknots(x,nknots,boundstab) ;

%	Create knots at sample quantiles.  If boundstab == 1, then nknot+2
%	knots are created and the first and last are deleted.  This
%	mitigates the extra variability of regression spline estimates near
%	the boundaries.
%	
%		INPUT (required)
%	x = independent variable.  (The knots are at sample quantiles of x.)
%	nknots = number of knots
%
%		INPUT (optional)
%	boundstab = parameter for boundary stability (DEFAULT is 0)
%
%	USAGE: knots = quantileknots(x,nknots,boundstab) ;
%
%
%	Last edit: 11/28/98
%
if nargin < 3 ;
	boundstab = 0 ;
end ;

n = length(x) ;
xsort = sort(x) ;   


loc = n*(1:nknots+2*boundstab)' ./ (nknots+1+2*boundstab) ;
knots=xsort(round(loc)) ;
knots=knots(1 + boundstab : nknots + boundstab) ;  
		%  REMOVE KNOTS NEAR BOUNDARIRES FOR
				%  STABILITY (= LOW VARIABILITY)
	
