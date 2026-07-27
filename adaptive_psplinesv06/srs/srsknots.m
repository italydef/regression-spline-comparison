function [yhat,beta,cp,imin,df,sigma2hat,xx,xy,xm ,...
		postvarbeta,postvaryhat,dfres] = ...
		srsknots(x,y,degree,knots,penwt,q);
%
%   function [yhat,beta,cp,imin,df,sigma2hat,xx,xy,xm, ...
%		postvarbeta,postvaryhat] = ...
%		srsknots(x,y,degree,knots,penwt,q);
%
%	Univariate smoother based on a penalized regression spline.
%	Uses a quadratic penalty.
%
%		INPUT
%	x = independent variable
%	y = dependent variable
%	degree = degree of the spline (default is 2) 
%	knots = set of knots (default is equally spaced quantiles)
%	penwt = candidates for the penalty weight - cp is used to choose
%		among these (default = logspace(-10,10,20))
%	q = penalty norm (default=2)
%
%		OUTPUT
%	yhat = estimated function
%	beta = regression spline coefficients
%	cp = cp at each value of penwt
%	imin = index of penwt that minimizes cp
%	df = effective degrees of freedom
%	sigma2hat = estimate of sigma squared
%	xx = xm'*xm where xm is the design matrix for the regression spline
%	xy = xm'*y
% 	xm - see above
%	postvarbeta = posterior variance of beta
%	postvaryhat = posterior variance of yhat
%
%
%
%	Last edit: Sept 19, 1997
%
if nargin < 3 ;
degree = 2 ;
end ;

if nargin < 4 ;
xsort = sort(x) ;                         %   Set knots if not input
loc = n*(1:nknots)' ./ (nknots+1) ;
knots=xsort(round(loc)) ;
xm=ones(n,1);
end ;

if nargin < 5 ;
penwt = logspace(-10,10,20) ;
end ;

if nargin < 6 ; 
q = 2; 
end ;

n = rows(x) ;                            
nknots = length(knots) ;
xm=ones(n,1);

	for i=1:degree ;
	xm = [xm x.^i] ;
	end ;
	
	for i=1:(nknots) ;
	xm = [xm ((x-knots(i)).^degree).*(x > knots(i))] ;
	end ;

	xx = xm'*xm;
	xy = xm'*y ;

m = max(size(penwt)) ;
beta = zeros(cols(xm),m) ;
yhat = zeros(n,m) ;
asr = zeros(m,1) ;
cp = asr ;
trsd = asr ;
if q==2;niter=ones(m,1);else;niter=[15;3*ones(m-1,1)];end;

	for i=1:m ;

 	id = diag([zeros(1,degree+1)  ones(1,nknots)]) ;
	
		for l=1:niter(i) ;
			if l == 1 ;
				if i==1;id2 = id ;else;
				id2 =diag([zeros(1,degree+1) ...
		                    lpwt(beta(degree+2:size(id,1),i-1)',q)]) ;
				end;
			else ;
			beta2 = beta(degree+2:size(id,1),i)' ;
			id2 = diag([zeros(1,degree+1) ...
		                    lpwt(beta2,q)]) ;
			end ;
		end ;


	b = inv(xx + penwt(i)*id) ;
	xxb = xx*b ;
	trsd(i) = trace(xxb) ;
	trsdsd(i) = trace(xxb*xxb) ;
	
	beta(:,i) = b * xy ;
	yhat(:,i) = xm*beta(:,i) ;
	asr(i) = mean((y-yhat(:,i)).^2) ;
		if i==1;;
		sigma2hat= n*asr(i)/(n-2*trsd(i)+trsdsd(i)) ;
		end ;
	cp(i) = asr(i) + 2*trsd(i)*sigma2hat/n ;
	cp(i) = n*asr(i)/sigma2hat + 2*trsd(i) - n ;
	
	end ;
imin = find(  (cp==min(cp))  ) ;
imin=min(imin) ;
b = inv(xx + penwt(imin)*id) ;

df= trsd ;
yhat = yhat(:,imin) ;
beta = beta(:,imin) ;

dfres = n-2*trsd(imin)+trsdsd(imin) ;
sigma2hat= n*asr(imin)/ dfres;

postvarbeta = sigma2hat*b ;
postvaryhat = diag(xm*postvarbeta*xm') ;

