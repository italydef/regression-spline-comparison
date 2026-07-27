function xm = addbasis(x,degree,nknots,z,xforknots,der,jder) ;

%	Outputs derivatives wrt jth x of the spline basis
%
%	USAGE:  xm = addbasis(x,degree,nknots,z,xforknots,der,jder) 
%
%		INPUT (REQUIRED)
%	x = matrix of independent variables (n by d, where d is the number
%		of independent variables
%
%		INPUT (OPTIONAL)
%	degree = degree of spline (default is 2)
%	nknots = number of knots for each independent variable (this is
%		d by 1.  If the input value is 1 by 1, then it is expanded
%		to a d by 1 constant vector.  (default is 10 for each
%		independent variable)
%	z = matrix of independent variable that enter linearly (can
%		be empty --- that is the default)
%	xforknots = x matrix used to create knots (can be the same as
%		x, but can be different --- default is xforknots = x)
%	der = order of the derivative (default = 0)
%	jder = derivative with respect to the jth component of x
%		(default is 1)
%
%
%		OUTPUT
%	xm = derivatives wrt jth x of the additive spline basis

%	Last edit:	1/6/99

[n,d] = size(x) ;

if nargin < 2 ;
	degree = 2 ;
end ;

if nargin < 3 ;
	nknots = 10 ;
end ;

if nargin < 4 ;
	z = [] ;
end ;

if nargin < 5 ;
	xforknots = x ;
end ;

if nargin < 6 ;
	der = 0 ;
end ;

if nargin < 7 ;
	jder = 1 ;
end ;

if der == 0 ;
	xm=ones(n,1);
	else ;
	xm = zeros(n,1) ;
end ;

if isempty(z) == 0 ;
	xm = [xm (der==0)*z] ;
end ;


	for j = 1:d ;
	knots = quantileknots(xforknots(:,j),nknots(j)) ; % Get knots at sample 
							% quantiles
	xmj = powerbasis(x(:,j),degree,knots,der) ;	% GET SPLINE BASIS
	if der > 0 & j ~= jder ;
		xmj = 0*xmj ;
	end ;
	xm = [xm xmj(:,2:size(xmj,2))] ;   
		%  xm is the "design matrix" of the regression spline
	end ;
