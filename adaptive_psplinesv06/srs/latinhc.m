function design = latinhc(n,k,lower,upper) ;
%	Returns a Latin hypercube sample
%
%		INPUT (required)
%	n = number of observations (rows)
%	k = number of variables (columns)
%
%		INPUT (optional)
%	lower = lower limit of each variable (default = 0)
%	upper = upper limit of each variable (default = 1)
%
%		OUTPUT
%	design = n by k matrix
%
%	The ith column of design has values in the interval [lower(i),upper(i)]
%
%	If lower and/or upper is scalar then it is expanded to a 
%	vector of constants
%
%	USAGE: design = latinhc(n,k,lower,upper) ;
%
%	Last edit:	July 3, 1998
%	
%	Example: design = latinhc(100,5,2,(3:7)) ;
%		This gives a 100 by 5 matrix.  The ith column has values
%		in the interval [2, 2+i]
%
if nargin < 4 ;
	upper = 1 ;
end ;
if nargin < 3 ;
	lower = 0 ;
end ;

if size(upper,2) > 1 ;
	upper = upper' ;
end ;
if size(lower,2) > 1 ;
	lower = lower' ;
end ;

if length(upper) == 1 ;
	upper = upper*ones(k,1) ;
end ;
if length(lower) == 1 ;
	lower = lower*ones(k,1) ;
end ;

range = upper-lower ;
design = rand(n,k)-1 ;
for i=1:k ;
	design(:,i) = (range(i)/n)*(design(:,i) + randperm(n)') + lower(i) ;
end ;
