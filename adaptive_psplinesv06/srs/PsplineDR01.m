function  fit = PsplineDR01(x,y,degree,nknots,penwt,boundstab);
%
%	Fits a P-spline to univariate x's with Demmler-Reinsch
%	algorithm.
%
%
%	USAGE:  fit = PsplineDR01(x,y,degree,nknots,penwt);
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
%		minimizing gcv (default is logspace(-12,12,100))
%	boundstab = parameter passed to quantileknots (see that program)
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
%	xm = design matrix of the spline (not returned if n > 1000)
%	xx = xm' * xm
%	a =  penwt(imin)
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
%	'llimitder',llimitder,'asr',asr) ;
%
%
%	CALLS: powerbasis01, quantileknots
%
%	Last edit: 7/26/99
%


n = size(x,1) ;

if (nargin < 6 | isempty(boundstab) == 1) ;
	boundstab = 0 ;
end ;

if (nargin < 5 | isempty(penwt) == 1) ;
penwt = logspace(-12,12,100) ;
end ;

if (nargin < 4 | isempty(nknots) == 1) ;
nknots = 20 ;
end ;

if (nargin < 3 | isempty(degree) == 1) ;
degree = 2 ;
end ;

stdx = std(x) ;
x = (x - mean(x)) ./ stdx ;
meany = mean(y) ;
y = y - meany ;

knots = quantileknots(x,nknots,boundstab) ;

n = length(x) ;

xm = powerbasis01(x,degree,knots) ;
xx = xm'*xm ;



id = [zeros(1,degree+1)  ones(1,nknots)] ;
D = diag(id) ;

%R=chol(xx + 10e-10*eye(size(xx,1)) ) ;


[R,p]=chol(xx) ;

if p > 0 ;
p
end ;



B = inv(R') ;
[U,C] = eig(B*D*B') ;
Z = xm*B'*U ;

Zy = Z'*y ;
ZZ=Z'*Z ;

m = length(penwt) ; 
beta = zeros(size(xm,2),m) ;
%yhat = zeros(n,m) ;
asr = zeros(m,1) ;
gcv = asr ;
trsd = asr ;
trsdsd = asr ;
dfres = asr ;
ssy = y'*y ;

	for i=1:m ;
	oneld = 1 ./ (1 + penwt(i)*diag(C)) ;
	trsd(i) = sum(oneld ) ;
	alpha(:,i) = (Zy) .* oneld ;
	asr(i) =  (ssy - 2*Zy'*alpha(:,i) + alpha(:,i)'*ZZ*alpha(:,i))/n ;
	trsdsd(i) = sum(oneld.*oneld) ;
	dfres(i) = n -2*trsd(i) + trsdsd(i) ;
	gcv(i) = asr(i) / (1-trsd(i)/n)^2 ;
	end ;

imin = min(find(  (gcv==min(gcv)) ) ) ;
a = penwt(imin) ;

dffit= trsd ;
alpha = alpha(:,imin) ;
beta = B'* U * alpha ;
dbeta = length(beta) ;
yhat = Z*alpha ;
sigma2hat = n*asr(imin) ./ dfres(imin) ;

oneld = 1 ./ (1 + penwt(imin)*diag(C)) ;
postvaralpha = sigma2hat*diag(oneld);
postvaryhat = (Z.*(Z*postvaralpha))*ones(dbeta,1) ;
postvarbeta = B'*U*postvaralpha*U'*B ;

Z = [] ;

xmder = xm ;


xmder(:,1) = 0*xmder(:,1) ;

for i = 2:degree+1 ;
	xmder(:,i) = (i-1)*(abs(xm(:,i))).^( (i-2)/(i-1)  ) .* (sign(x)).^(i-2);
end ;

for i=degree+2:degree+1+nknots ;
	xmder(:,i) = degree* (abs(xm(:,i))).^(  (degree-1)/degree ) ...
		.* (xm(:,i) > 0);
end ;



yhatder = xmder*beta ;
postvaryhatder = ((xmder*postvarbeta).*xmder)*ones(dbeta,1) ;

clear xmder ;

yhat = yhat + meany ;
ulimit = yhat + 2*sqrt(postvaryhat) ;
llimit = yhat - 2*sqrt(postvaryhat) ;
ulimitder = yhatder + 2*sqrt(postvaryhatder) ;
llimitder = yhatder - 2*sqrt(postvaryhatder) ;


if n < 1001 ;
fit = struct('yhat',yhat,'beta',beta,'gcv',gcv,'imin',imin,'dffit',dffit, ...
	'knots',knots,'postvarbeta',postvarbeta,'postvaryhat',postvaryhat, ...
	'xm',xm,'xx',xx,'a',a,'penwt',penwt,'sigma2hat',sigma2hat, ...
	'dfres',dfres, ...
	'yhatder',yhatder,'postvaryhatder',postvaryhatder, ...
	'ulimit',ulimit,'llimit',llimit,'ulimitder',ulimitder, ...
	'llimitder',llimitder, ...
	'asr',asr) ;
else ;
fit = struct('yhat',yhat,'beta',beta,'gcv',gcv,'imin',imin,'dffit',dffit, ...
	'knots',knots,'postvarbeta',postvarbeta,'postvaryhat',postvaryhat, ...
	'xx',xx,'a',a,'penwt',penwt,'sigma2hat',sigma2hat, ...
	'dfres',dfres, ...
	'yhatder',yhatder,'postvaryhatder',postvaryhatder, ...
	'ulimit',ulimit,'llimit',llimit,'ulimitder',ulimitder, ...
	'llimitder',llimitder, ...
	'asr',asr) ;

end ;
