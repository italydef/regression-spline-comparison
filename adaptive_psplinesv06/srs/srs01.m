function  fit = srs01(x,y,degree,nknots,penwt,boundstab);
%
%	Fits a regression spline to univariate x's with a quadratic penalty
%
%	USAGE:  fit = srs01(x,y,degree,nknots,penwt);
%
%
%		INPUT - REQUIRED
%	x = independent variable (univariate)
%	y = response (same length as x)
%
%		INPUT - OPTIONAL
%	degree = degree of the spline (default is 2)
%	nknots = number of knots (default 20)
%	penwt = trial values of the penalty weight (one is chosen by 
%		minimizing gcv (default is logspace(-10,5,31))
%	
%		OUTPUT
%	Returns a structure "fit" with the following components:
%
%	yhat = fitted values (same length as x)
%	beta = spline coefficients
%	gcv = vector of gcv values (same length as penwt)
%	imin = index of penwt where minimum of gcv occurs
%	dffit = degrees of freedom of the smoother (same length as penwt)
%	knots = knots (length = nknots)
%	postvarbeta = posterior variance of beta (from a Bayesian model)
%	postvaryhat = posterior variance of yhat (from a Bayesian model)
%	xm = design matrix of the spline
%	xx = xm' * xm
%	a = alpha = penwt(imin)
%	penwt
%	sigma2hat
%	dfres
%	varbeta
%	stdyhat
%	yhatder
%	postvaryhatder
%	ulimit
%	llimit
%	ulimitder
%	llimitder
%
%		That is, 
%	fit = struct( ...
%	'yhat',yhat,'beta',beta,'gcv',gcv,'imin',imin,'dffit',dffit, ...
%	'knots',knots,'postvarbeta',postvarbeta,'postvaryhat',postvaryhat, ...
%	'xm',xm,'xx',xx,'a',a,'penwt',penwt,'sigma2hat',sigma2hat, ...
%	'dfres',dfres,'varbeta',varbeta,'stdyhat',stdyhat, ...
%	'yhatder',yhatder,'postvaryhatder',postvaryhatder, ...
%	'ulimit',ulimit,'llimit',llimit,'ulimitder',ulimitder, ...
%	'llimitder',llimitder) ;
%
%
%	CALLS: powerbasis, quantileknots
%
%	Last edit: 7/16/99
%
n = size(x,1) ;

if (nargin < 6 | isempty(boundstab) == 1) ;
	boundstab = 0 ;
end ;

if (nargin < 5 | isempty(penwt) == 1) ;
penwt = logspace(-10,10,30) ;
end ;

if (nargin < 4 | isempty(nknots) == 1) ;
nknots = 20 ;
end ;

if (nargin < 3 | isempty(degree) == 1) ;
degree = 2 ;
end ;

stdx = std(x) ;
x = (x - mean(x)) ./ stdx;

knots = quantileknots(x,nknots,boundstab) ;
n = length(x) ;
xm=ones(n,1);
xm = powerbasis01(x,degree,knots) ;
xx = xm'*xm ;
xy = xm'*y ;
id = diag([zeros(1,degree+1)  ones(1,nknots)]) ;
m = length(penwt) ; 
beta = zeros(size(xm,2),m) ;
yhat = zeros(n,m) ;
asr = zeros(m,1) ;
gcv = asr ;
trsd = asr ;
trsdsd = asr ;
dfres = asr ;
ssy = y'*y ;

	for i=1:m ;
	xxp = xx + penwt(i)*id ;
	xxb = xx/xxp ;
	trsd(i) = trace(xxb) ;
	beta(:,i) = xxp\xy ;
	asr(i) =  (ssy - 2*xy'*beta(:,i) + beta(:,i)'*xx*beta(:,i))/n ;
	trsdsd(i) = trace(xxb*xxb) ;
	dfres(i) = n -2*trsd(i) + trsdsd(i) ;
	gcv(i) = asr(i) / (1-trsd(i)/n)^2 ;

	end ;

imin = min(find(  (gcv==min(gcv)) ) ) ;

dffit= trsd ;
sigma2hat = n*asr(imin) ./ dfres(imin) ;
ixxp = inv(xx + penwt(imin)*id) ;
beta = beta(:,min(imin)) ;
yhat = xm*beta ;
postvarbeta = sigma2hat*ixxp ;
xxp = xx + penwt(imin)*id ;
postvaryhat = sigma2hat * (xm.*(xm/xxp))*ones(length(beta),1) ;
varbeta = sigma2hat * ixxp * xx * ixxp ;
stdyhat = sqrt(   (xm.*(xm*varbeta))*ones(length(beta),1)   ) ;
a = penwt(imin) ;


xmder = xm ;
xmder(:,1) = 0*xmder(:,1) ;
for i = 2:degree+1 ;
xmder(:,i) = (i-1)*(abs(xm(:,i))).^( (i-2)/(i-1)  ) .* (sign(x)).^(i-2);
end ;

for i=degree+2:degree+1+nknots ;
xmder(:,i) = degree* (abs(xm(:,i))).^(  (degree-1)/degree ) ...
	.* (xm(:,i) > 0);
end ;

yhatder = xmder*beta / stdx ;
postvaryhatder = (sigma2hat/ (stdx^2)) * ... 	
	(xmder.*( xmder/xxp ))  *ones(length(beta),1) ;


ulimit = yhat + 2*sqrt(postvaryhat) ;
llimit = yhat - 2*sqrt(postvaryhat) ;


ulimitder = yhatder + 2*sqrt(postvaryhatder) ;
llimitder = yhatder - 2*sqrt(postvaryhatder) ;



fit = struct('yhat',yhat,'beta',beta,'gcv',gcv,'imin',imin,'dffit',dffit, ...
	'knots',knots,'postvarbeta',postvarbeta,'postvaryhat',postvaryhat, ...
	'xm',xm,'xx',xx,'a',a,'penwt',penwt,'sigma2hat',sigma2hat, ...
	'dfres',dfres,'varbeta',varbeta,'stdyhat',stdyhat, ...
	'yhatder',yhatder,'postvaryhatder',postvaryhatder, ...
	'ulimit',ulimit,'llimit',llimit,'ulimitder',ulimitder, ...
	'llimitder',llimitder, ...
	'asr',asr) ;
