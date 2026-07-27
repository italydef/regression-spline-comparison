function df = diffspline(x,knots,beta,degree,order) ;
%
%	Differentiates a regression spline with given knots, degree,
%	and beta coefficient vector.
%
%		INPUT (required)
%	x = points where the derivative will be evaluated
%	knots = knots of the spline
%	beta = coefficient vector using polynomial and truncated power 
%		basis
%	degree = degree of the spline
%	order = order of the derivative
%
%		OUTPUT
%	df = derivative of the spline
%
%	Last edit: May 4, 1998
%
%	Copyright: David Ruppert
%

coeff=(fact(degree)/fact(degree-order))*ones(length(beta),1) ;

for i=1:degree;
	if i < order+1 ;
		coeff(i) = 0 ;
	elseif i <= degree ;
		coeff(i) = fact(i-1)/fact(i-order-1) ;
	end ;
end ;

beta = coeff.*beta ;
betader = beta(order+1:length(beta)) ;
xm = powerbasis(x,degree-order,knots) ;
df = xm*betader ;
