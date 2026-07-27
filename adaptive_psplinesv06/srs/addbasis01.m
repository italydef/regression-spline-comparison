function basis = addbasis01(x,degree,nknots,z,xforknots,der,jder,nsubknots) ;
%
%	Copied from addbasis.m
%
%	Outputs derivatives wrt jth x of the spline basis
%
%	USAGE:  xm = addbasis01(x,degree,nknots,z,xforknots,der,jder,nsubknots) 
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
%	nsubknots = number of subknots for each x variable (scalar)
%
%
%		OUTPUT
%	xm = derivatives wrt jth x of the additive spline basis

%	Last edit:	7/12/99

[n,d] = size(x) ;

if (nargin < 2 | isempty(degree) == 1 ) ;
	degree = 2 ;
end ;

if ( nargin < 3 | isempty(nknots) == 1 ) ;
	nknots = 10 ;
end ;

if nargin < 4  ;
	z = [] ;
end ;

if ( nargin < 5 | isempty(xforknots) == 1 );
	xforknots = x ;
end ;

if ( nargin < 6 | isempty(der) == 1 ) ;
	der = 0 ;
end ;

if ( nargin < 7 | isempty(jder) == 1 );
	jder = 1 ;
end ;

if ( nargin < 8 | isempty(nsubknots) == 1 ) ;
	nsubknots = 0 ;
end ;

if der == 0 ;
	xm=ones(n,1);
	else ;
	xm = zeros(n,1) ;
end ;

if isempty(z) == 0 ;
	xm = [xm (der==0)*z] ;
end ;

if length(nknots) == 1 ;
	nknots = ones(d,1)*nknots ;
end ;
maxnknots = max(nknots) ;
knots = zeros(maxnknots,d) ;


if nsubknots > 0 ;
	subknots = zeros(nsubknots,d) ;
else ;
	subknots = [] ;
end ;


for j = 1:d ;
	xunique = unique(xforknots(:,j)) ;

	knots(1:nknots(j),j) = quantileknots(xunique, ...
		nknots(j)) ;

	if nsubknots > 0 ;
		xunique2 = xunique(  (xunique>knots(1,j)) & ...
			(xunique<knots(nknots(j),j)) ) ;
	
		subknots2 = quantileknots(xunique2,nsubknots-2) ;

		subknots(:,j) = [knots(1,j); subknots2; knots(nknots(j),j)] ;
	end ;


	xmj = powerbasis(x(:,j),degree,knots(1:nknots(j),j), ...
			der) ;	% GET SPLINE BASIS
	if der > 0 & j ~= jder ; 
		xmj = 0*xmj ;
	end ;
	xm = [xm xmj(:,2:size(xmj,2))] ;   
		%  xm is the "design matrix" of the regression spline
end ;

basis = struct('xm',xm,'knots',knots,'subknots',subknots) ;
