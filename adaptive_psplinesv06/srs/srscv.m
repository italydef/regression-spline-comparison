function [yhat,beta,gcv,imingcv,df,knots,postvarbeta,postvaryhat,xm,xx,a, ...
	penwt,sigma2hat,cv,imincv,asr] = srs(x,y,degree,nknots,penwt);
%
%	Fits a regression spline to univariate x's with a quadratic penalty
%
%	Similar to srs but also computes cv.  Is slow and should only
%	be used with small data sets.
%
%	USAGE: [yhat,beta,gcv,imin,df,knots,postvarbeta,postvaryhat,xm,xx,a, ...
%	penwt,sigma2hat] = srs(x,y,degree,nknots,penwt);
%
%
%		INPUT - REQUIRED
%	x = independent variable (univariate)
%	y = response (same length as x)
%
%		INPUT - OPTIONAL
%	degree = degree of the spline (default is 2)
%	nknots = number of knots (default is min(n/3,40))
%	penwt = trial values of the penalty weight (one is chosen by 
%		minimizing gcv (default is logspace(-10,5,31))
%	
%		OUTPUT
%	yhat = fitted values (same length as x)
%	beta = spline coefficients
%	imin = index of penwt where minimum of gcv occurs
%	cv =  vector of cv values (same length as penwt)
%	gcv = vector of gcv values (same length as penwt)
%	df = degrees of freedom of the smoother (same length as penwt)
%	knots = knots (length = nknots)
%	postvarbeta = posterior variance of beta (from a Bayesian model)
%	postvaryhat = posterior variance of yhat (from a Bayesian model)
%	xm = design matrix of the spline
%	xx = xm' * xm
%	a = alpha = penwt(imin)
%
%	CALLS: powerbasis, quantileknots
%
%	Last edit: 1/23/98
%
n = size(x,1) ;
if nargin < 5 ;
penwt = logspace(-6,7,30) ;
end ;

if nargin < 4 ;
nknots = ceil(min(n/3,40)) ;
end ;

if nargin < 3 ;
degree = 2 ;
end ;

knots = quantileknots(x,nknots) ;
n = length(x) ;
xm=ones(n,1);
xm = powerbasis(x,degree,knots) ;
xx = xm'*xm ;
xy = xm'*y ;
id = diag([zeros(1,degree+1)  ones(1,nknots)]) ;
m = length(penwt) ; 
beta = zeros(size(xm,2),m) ;
yhat = zeros(n,m) ;
asr = zeros(m,1) ;
cv = asr ;
gcv = asr ;
trsd = asr ;
ssy = y'*y ;

	for i=1:m ;
	b = inv(xx + penwt(i)*id) ;
	xxb = xx*b ;
	trsd(i) = trace(xxb) ;
	beta(:,i) = b * xy ;
	asr(i) =  (ssy - 2*xy'*beta(:,i) + beta(:,i)'*xx*beta(:,i))/n;
		if i==1;
		trsdsd = trace(xxb*xxb) ;
		sigma2hat= n*asr(i)/(n-2*trsd(i)+trsdsd) ;
		end ;
	gcv(i) = asr(i) / (1-trsd(i)/n)^2 ;
	cv(i) = sum((  (y-xm*beta(:,i)) ./ (1-diag(xm*b*xm')) ).^2 )./n ;

	end ;

imingcv = min(find(  (gcv==min(gcv)) ) ) ;
imincv = min(find(  (cv==min(cv)) ) ) ;

df= trsd ;
sigma2hat = n*asr(imincv) ./ (n - df(imincv)) ;
b = inv(xx + penwt(imincv)*id) ;
beta = beta(:,min(imincv)) ;
yhat = xm*beta ;
postvarbeta = sigma2hat*b ;
postvaryhat = sigma2hat*(xm.*(xm*b))*ones(length(beta),1) ;
a = penwt(imincv) ;
