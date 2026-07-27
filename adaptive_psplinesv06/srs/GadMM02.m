function results = GadMM01(x1,y,x2,degree,nknots,niter,sigfactor) ;
%
%	Fits a generalized additive mixed model by MCMC
%
%	CURRENTLY ONLY FITS LOGISTIC REGRESSION
%
%		INPUT (required)
%	y = response (modeled as conditionally Binomial or
%		Poisson (n by 1 dimensional - POISSON NOT YET IMPLEMENTED)
%	x1 = matrix of covariates modeled additively (n by m1 dimensional)
%
%		INPUT (optional)
%	niter = number of iterations of the Markov Chain sampler
%		(default = 1000)
%	x2 = matrix of covariate modeled linearly (default = [])
%
%	Last edit: 6/16/2000
%

	%	Fill in missing arguments with defaults

if (nargin < 7 | isempty(sigfactor) == 1) ;	
	sigfactor = .25 ;
end ;
if (nargin < 6 | isempty(niter) == 1) ;	
	niter = 1000 ;
end ;

if (nargin < 5 | isempty(nknots) == 1);
	nknots = 5 ;
end ;

if (nargin < 4 | isempty(degree) == 1) ;
	degree = 1 ;
end ;

if nargin < 3  ;
	x2 = [] ;
end ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%		Set up sampler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=length(y) ;
m1=size(x1,2) ;		%	m1 = number of covariates modeled additively
m2=size(x2,2) ;		%	m2 = number of covariates modeled linearly
fitprel = LogSpl01(x1,y,degree,nknots,x2) ;	%	Get preliminary fit
basis = addbasis01(x1,degree,nknots,x2) ;
xm = basis.xm ;

nknots = fitprel.nknots ;
stdbeta = sqrtm(fitprel.varmatrix) ;
sigcand = sigfactor*stdbeta  ;
beta = fitprel.beta ;
gcv=fitprel.gcv(fitprel.imin) ;
sigu = ones(m1,1)/gcv ;		%	Preliminary estimate of prior variances

m = size(xm,2) ;	%	m = dimension of beta (beta = length of ALL
			%		regression parameters


A = .01 ;		%%%%%	Define hyperparameters
B = .01 ;
taudiff = 10^(-8) ;
tauvect = ones(1+m2,1)*taudiff ; 
tauu = 0*nknots ;

for i=1:m1 ;
		pointer = 2+m2 + sum(nknots(1:i-1)) + (i-1)*degree ;
		coeff = beta( pointer : pointer - 1 + nknots(i) ) ;

		tauu(i) =  gamrnd(A+nknots(i)/2, 1/(B + sum(coeff.^2)/2) ) ;
	tauvect = [tauvect;ones(degree,1)*taudiff;ones(nknots(i),1)*tauu(i) ] ;
end;


sigu = 1./tauu ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%		Set up storage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
betastore = zeros(niter,m) ;
sigustore = zeros(niter,m1) ;
acceptbeta = zeros(niter-1,1) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%		Start of Markov Chain Sampler 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for iter=1:niter-1 ;
betastore(iter,:) = beta' ;
sigustore(iter,:) = sigu' ;


		%	Metropolis Hastings step for beta 

	
	unif = rand(1,1) ;
	candbeta = beta + sigcand*randn(m,1) ;
	logquotbeta = y'*xm*candbeta - sum( log(1+exp(xm*candbeta)) ) ...
		- tauvect.*(candbeta.^2)/2 ...
		 - (y'*xm*beta - sum( log(1+exp(xm*beta)) ) - ...
		tauvect.*(beta.^2)/2) ;


	indbeta = (log(unif) < sum(logquotbeta)) ;
	beta = indbeta.*candbeta + (1-indbeta).* beta ;
	acceptbeta(iter,:) = indbeta' ;


	 	% Generate sigu as inverse gammas
	tauvect = ones(1+m2,1)*taudiff ;

	for i=1:m1 ;
		pointer = 2+m2 + sum(nknots(1:i-1)) + (i-1)*degree ;
		coeff = beta( pointer : pointer - 1 + nknots(i) ) ;

		tauu(i) =  gamrnd(A+nknots(i)/2, 1/(B + sum(coeff.^2)/2) ) ;
		sigu(i) = 1/tauu(i) ;
		tauvect = [tauvect;ones(degree,1)*taudiff; ...
		ones(nknots(i),1)*tauu(i) ] ;
	end;	%	End "for i=1:m1"



end ;		%	End "for iter=1:niter-1"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%		End of Markov Chain Sampler 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


betastore(niter,:) = beta' ;
sigustore(niter,:) = sigu' ;
ngrid = 200 ;
xbar = ones(ngrid,1)*mean(x1) ;
if isempty(x2) == 0 ;
	zbar = ones(ngrid,1)*mean(x2) ;
else ;
	zbar = [] ;
end ;

xgrid = zeros(ngrid,m1) ;
start = min([1000 floor(.1*niter)]) ;	%	Define burn-in period

meanbeta = mean(betastore(start:niter,:))' ;

for j=1:m1 ;
	xj = xbar ;
	xderj = 0*xbar ;
	xj(:,j) = linspace(min(x1(:,j)),max(x1(:,j)),ngrid)' ;
	xderj(:,j) = xj(:,j) ;
	xgrid(:,j) = xj(:,j) ;


	basis = addbasis01(xj,degree,nknots,zbar,x1,[],0) ;
	xmj = basis.xm ;
	mhat(:,j)=xmj*meanbeta;
	basis = addbasis01(xderj,degree,nknots,zbar,x1,1,j) ;	
	xmderj = basis.xm ;
	mhatder(:,j) = xmderj*meanbeta  ;
end ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%		Put results into a structure 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = struct('beta',betastore,'sigu',sigustore, ...
	'meanbeta',meanbeta,'mhat',mhat,'mhatder', mhatder, ...
	'acceptbeta',mean(acceptbeta),'fitprel',fitprel,'xgrid',xgrid, ...
	'xm',xm) ;
