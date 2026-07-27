function [yhat,beta,gcv,iminnormp,df,diagid,a,normp]=...
	srslpauto(x,y,degree,nknots,penwt);
%
%   Fits a regression spline to univariate x's with an L_p norm penalty
%
%
%
%   Last edited:  1/15/98
%
%		INPUT - REQUIRED
%	x = independent variable (univariate)
%	y = response (same length as x)
%
%		INPUT - OPTIONAL
%	degree = degree of the spline (default is 2)
%	nknots = number of knots (default is 20)
%	normp = p in the norm, e.g., normp = 2 for L_2 norm (default = 2)
%	penwt = vector of trial values of penalty weight
%		(default=logspace(-5,10,30)' 
%
%		OUTPUT
%	yhat = fitted values (same length as x)
%	beta = regression spline coefficients (length = 1 + degree + nknots)
%	gcv = gcv values (one for each value of penwt)
%	imin = index that minimizes gcv
%	df = vector of degrees of freedom (same length as penwt)
%	diagid = diag of the effective weighting matrix
%	a = alpha = penwt(imin)
%
%
%	CALLS:	powerbasis, lpwt
%
%	Copyright: David Ruppert

if nargin < 3 ;
degree = 2 ;
end ;

if nargin < 4 ;
nknots = 20 ;
end ;

if nargin < 5 ;
penwt = logspace(-11,10,40)' ;
end ;

n = length(x) ;
m = length(penwt) ;

normp =[.5 ;1;1.5;2];
np = length(normp) ;
gcv = zeros(m,np) ;
yhat = zeros(n,np) ;
beta = zeros(1+degree+nknots,np) ;
imin = gcv ;
df = zeros(m,np) ;

for i = 1:np;
[yhat(:,i),beta(:,i),gcv(:,i),imin(i),df(:,i),diagid,a] = ...
		srslp(x,y,degree,normp(i),nknots,penwt);
gcvmin(i) = gcv(imin(i),i) ;
end ;
gcvmin
iminnormp = min(find(gcvmin == min(gcvmin))) ;
normp = normp(iminnormp) ;
yhat = yhat(:,iminnormp) ;
beta = beta(:,iminnormp) ;

