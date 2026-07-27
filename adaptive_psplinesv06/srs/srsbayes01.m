function  fit = srsbayes01(x,y,degree,nknots,niter);
%
%	Gibbs sampler Bayesian analysis of a univariate regression
%	spline
%
%	USAGE:  fit = srsbayes01(x,y,degree,nknots,niter,boundstab);
%
%
%		INPUT - REQUIRED
%	x = independent variable (univariate)
%	y = response (same length as x)
%
%		INPUT - OPTIONAL
%	degree = degree of the spline (default is 2)
%	nknots = number of knots (default 20)
%	niter  = number of iterations of the Gibbs sampler (default is 2000)
%	boundstab = parameter passed to quantileknots (see that program)
%			(default = 0)
%	
%		OUTPUT
%	Returns a structure "fit" with the following components:
%
% 	fit = struct('beta',beta,'xmatrix',xm,'sigma2e',sigma2e,
%	'sigma2b',sigma2b,'knots',knots,'xgrid',xgrid,
%	'yhatgrid',yhatgrid, 'meanyhatgrid',meanyhatgrid, 
%	'varyhatgrid',varyhatgrid) ;

%	beta	= niter realizations of beta
%	xmatrix = xm = design matrix for the spline
%	sigma2e = niter realizations of sigma2e
%	sigma2b = niter realizations of sigma2b
%	knots 	= knots used in the spline (these are fixed and do not
%			vary between realizations)
%	xgrid	= grid of 100 points from min(x) to max(x)
%	yhatgrid = niter realization of yhat on xgrid
%	meanyhatgrid	= mean of yhatgrid
%	varyhatgrid	= variance of yhatgrid
%	varmeanyhatgrid = estimate of var(E(yhat|var components))
%	meanvaryhatgrid = estimate of E(var(yhat| var components))
%
%
%	Edit: 3/6/01
%

n = length(x) ;

if (nargin < 5 | isempty(niter) == 1) ;
	niter = 2000 ;
end ;

if (nargin < 4 | isempty(nknots) == 1) ;
nknots = 20 ;
end ;

if (nargin < 3 | isempty(degree) == 1) ;
	degree = 2 ;
end ;

stdx = std(x) ;
meanx = mean(x) ;
x = (x - meanx) ./ stdx ;
meany = mean(y) ;
y = y - meany ;

knots = quantileknots(x,nknots,0) ;


xm = powerbasis01(x,degree,knots) ;
xx = xm'*xm ;

id = [zeros(1,degree+1)  ones(1,nknots)] ;
D = diag(id) ;

R=chol(xx + 1000*eps*eye(size(xx,1)) ) ;
B = inv(R') ;
[U,C] = eig(B*D*B') ;
Z = xm*B'*U ;

Zy = Z'*y ;
ZZ=Z'*Z ;
fitprel = PsplineDR02(x,y,degree,nknots);
lambda =fitprel.penwt(fitprel.imin) ;
beta = fitprel.beta ;
sigma2hat = fitprel.sigma2hat ;

beta=zeros(1+degree+nknots,niter) ;
alphamean = beta ;
alpha = beta ;
varbeta = beta ;
sigma2e = zeros(1,niter) ;
sigma2b = sigma2e ;
varyhatgrid = zeros(100,niter) ;

Ab = .01 ;
Bb = .01 ;
Ae = .01;
Be = .01 ;
ssy = y'*y ;

xgrid=linspace(min(x),max(x),100)' ;
xgridm=powerbasis01(xgrid,degree,knots) ;

for i=1:niter;

oneld = 1 ./ (1 + lambda*diag(C)) ;
alphamean(:,i) = (Zy) .* oneld ;

postSDalpha = sqrt(sigma2hat)*diag(sqrt(oneld)) ;
postvaralpha = (sigma2hat)*diag((oneld)) ;


alpha(:,i) = alphamean(:,i)  + postSDalpha*randn(1+degree+nknots,1)  ;
beta(:,i) = B'* U * alpha(:,i) ;
varyhatgrid(:,i) = diag(xgridm*B'*U*postvaralpha*U'*B*xgridm') ;

ssr(i) =  ssy - 2*Zy'*alpha(:,i) + alpha(:,i)'*ZZ*alpha(:,i) ;
b = beta(degree+2:1+degree+nknots,i) ;
ssb = sum(b.^2) ;
gama= gamrnd(Ab+nknots/2, 1/(Bb + ssb/2) ) ;
sigma2b(i) = 1/gama ;
gama= gamrnd(Ae+n/2, 1/(Be + ssr(i)/2) ) ;
sigma2e(i) = 1/gama ;
lambda = sigma2e/sigma2b ;
beta(1,i) = beta(1,i)+meany ;

sigma2hat = sigma2e(i) ;

end ;

alphameanmean = mean(alphamean')' ;

betamean = B'*U*alphamean ;

betamean(1,:) = betamean(1,:) + meany*ones(1,niter) ;

yhatgrid = xgridm*beta  ;
yhatmeangrid = xgridm*betamean ;
yhatmeanmeangrid = mean(yhatmeangrid')' ;

varmeanyhatgrid  =var(yhatmeangrid')' ;

meanvaryhatgrid =mean(varyhatgrid')';

varyhatgrid = varmeanyhatgrid + meanvaryhatgrid ;

xgrid = meanx + stdx*xgrid ;

fit = struct('beta',beta,'xmatrix',xm,'sigma2e',sigma2e,'sigma2b',sigma2b,...
'knots',knots,'xgrid',xgrid,'yhatgrid',yhatgrid, ...
'meanyhatgrid',yhatmeanmeangrid, ...
'varyhatgrid',varyhatgrid,'varmeanyhatgrid',varmeanyhatgrid, ...
'meanvaryhatgrid',meanvaryhatgrid) ;



