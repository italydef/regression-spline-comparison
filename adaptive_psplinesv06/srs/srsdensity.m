function [yhat,x,beta,gcv,imin,df] = srsdensity(y,min_y,max_y,degree,nknots,penwt);
%
%   Regression spline density estimate
%
%   USAGE:  [yhat,x,beta,gcv,imin,df] = srs(y,min_y,max_y,degree,nknots,penwt)
%
%		OUTPUT
%	yhat	= estimated density at values of x
%	x	= estimation points 
%			(points where the density has been estimated)
%	gcv	= gcv values at values of penwt
%	imin	= index of penwt that minimizes gcv
%	df	= degree of freedom as value of penwt
%
%		INPUT (required)
%	y	= sample
%
%		INPUT (optional)
%	min_y	= smallest value at which the density is estimated
%	max_y	= largest value at which the density is estimated
%	degree 	= degree of the spline
%	nknots	= number of knots
%	penwt	= penalty weight
%
if nargin < 2 ;
min_y = min(y) - range(y)/3 ;
end ;

if nargin < 3 ;
max_y = max(y) + range(y)/3 ;
end ;

if nargin < 4 ;
degree = 2 ;
end ;

if nargin < 5 ;
nknots = 15 ;
end ;

if nargin < 6 ;
penwt = logspace(-10,5,35)' ;
end ;

nbin = 200 ;
n = length(y);
x=linspace(min_y,max_y,nbin)';
y = hist(y,x)' ;
supp_y = max_y - min_y;
coeff = nbin/(n*supp_y) ;
y = y .* coeff ;
knots=linspace(min_y,max_y,nknots+2)' ;
knots = knots(2:nknots+1) ;
xm = powerbasis(x,degree,knots) ;

id = diag([zeros(1,degree+1)  ones(1,nknots)]) ;
m = length(penwt) ;
dim_beta = size(xm,2) ;
beta = zeros(dim_beta,m) ;
yhat = zeros(nbin,m) ;
asr = zeros(m,1) ;
gcv = asr ;
trsd = asr ;

	for i=1:m ;
		if i==1 ;
		start = -log(max_y-min_y) ;
		beta0=[start;zeros(dim_beta-1,1)] ;
		niter = 200 ;
		else ;
		beta0 = beta(:,i-1) ;
		niter = 6 ;
		end ;

		for iter=1:niter ;
		mu = exp(xm*beta0) ;
		wxm = ((1./mu)*ones(1,size(xm,2))).*xm ;
		b = inv(xm'*wxm + penwt(i)*id) ;
		beta0 = b*(wxm'*(y-mu) + xm'*wxm*beta0) ;
		end ;
	xxb = wxm'*xm*b ;
	trsd(i) = trace(xxb) ;
	
	
	beta(:,i) = beta0 ;
	yhat(:,i) = exp(xm*beta(:,i)) ;
	asr(i) =  norm( (y-yhat(:,i)) )^2 ;

	end ;
gcv = asr./ (1-trsd./nbin).^2 ;

imin = min(find(  (gcv==min(gcv)) ) ) ;
df= trsd ;
yhat = yhat(:,imin) ;
beta = beta(:,imin) ;
