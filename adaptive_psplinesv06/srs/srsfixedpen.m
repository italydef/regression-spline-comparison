function [gcv,yhat,beta,asr,df,postvarbeta,dfres] = srsfixedpen(alpha,x,y,knots,subknots,...
	xx,xy,xm,sigma2hat,degree);
%
%   Fits a regression spline to univariate x's with a quadratic penalty
%   The penalty weight is fixed but may depend on the observation, i.e.,
%	may be local.
%
%   USAGE:  [gcv,yhat,beta,asr,df,postvarbeta] = 
%		srsfixedpen(alpha,x,y,knots,subknots,...
%			xx,xy,xm,sigma2hat,degree);
%
%
%   Defaults
%
%	Last edit: 7/8/99
%
if nargin < 10 ;
degree = 2 ;
end ;

n = rows(x) ;
nknots=length(knots) ;
penalty = interp1(subknots,alpha,knots,'linear') ;
binv = xx + diag([zeros(degree+1,1);exp(penalty)]) ;

b = inv(binv) ;
xxb = binv\xx ;
trsd = trace(xxb) ;
beta = binv\xy ;

yhat = xm*beta ;
asr = mean((y-yhat).^2) ;
trsdsd = trace(xxb*xxb) ;
df= trsd ;
dfres = n -2*trsd + trsdsd ;

gcv = asr / (1-trsd/n)^2 ;


postvarbeta = sigma2hat*b ;
