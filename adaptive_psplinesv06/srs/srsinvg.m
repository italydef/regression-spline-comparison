function [yhat,beta,logpost,imax,df,sigma2hat,gcv,pwt,s1,asr,trdf,yhatgcv] = ...
		 srsinvg(x,y,nknots,degree,q,pwt,r,r1);
%
%   Fits a regression spline using a Bayesian approach.  Uses the
%   Box-Tiao family for the prior on the regression coefficients
%   and inverse gamma priors on the sigmas.  Finds MAP estimate.
%   No GIBBS sampling.
%
%	INPUTS:
%	x = independent variable
%	y = dependent variable
%	nknots = number of knots (default = 40)
%	degree = degree of spline (default = 1)
%	q = norm parameter for penalty ( default = 2)
%	pwt = vector of candidate penalty weights (default=logspace(-10,10,21)')
%
%	OUTPUT:
%	yhat = spline estimation at each x
%	beta = coefficients of the spline model
%	logpost = log posterior as a function of pwt (rows) and s1 (columns)
%	imax,jmax = indices maximizing logpost
%	df = degrees of freedom
%	sigma2hat = estimate of sigma
%	pwt = same as input or default value if not input
%	s1 =
%	
%
%   	Last edited:  7/9/97
%
%
%
%     CALLS: lpwt
%
%

if nargin < 8 ;
r1 = 1/2 ;
end ;
 
if nargin < 7 ;
r = 1/2 ;
end ;

if nargin < 6 ;
pwt = logspace(-10,10,21)' ;
end ;

if nargin < 5 ;
q = 2 ;
end ;

if nargin < 4 ;
degree = 1 ;
end ;

if nargin < 3 ;
nknots = 40 ;
end ;

n = size(x,1) ;


xsort = sort(x) ;                         %   Set knots
loc = n*(1:nknots)' ./ (nknots+1) ;
knots=xsort(round(loc)) ;

[yhat,betaprel,cp,imin,df,sigma2hat0,xx,xy,xm] = ...
		srsknots(x,y,degree,knots,pwt);  % Initial smooth to
						%    estimate sigma
q3beta = median(abs(betaprel)) ;

s=sigma2hat0 ;
 
s1= q3beta^2 ;

if q == .5 ;
s1 = s1/(11.25)^2 ;
elseif q == 1 ;
s1 = s1/(1.38)^2 ;
elseif q==2 ;
s1 = s1/(.674)^2 ;
end ;

w = ones(n,1) ;
xm=ones(n,1);

	for i=1:degree ;
	xm = [xm x.^i] ;
	end ;
	
	for i=1:(nknots) ;
	xm = [xm ((x-knots(i)).^degree).*(x > knots(i))] ;
	end ;
	
xmw = xm.*(w*ones(1,size(xm,2))) ;
yw = y.*w ;
xx = xmw'*xmw ;
xy = xmw'*yw ;

id = diag([zeros(1,degree+1)  ones(1,nknots)]) ;
m = max(size(pwt)) ;         %  pwt is the vector of penalty weights
beta = zeros(size(xm,2),m) ;   %  Set vectors for storage
yhat = zeros(n,m) ;          %  Fitted values as a function of pwt
asr = zeros(m,1) ;           %  Average squared residuals
logpost = zeros(m,length(s1)) ;              %  log posterior
logposta = logpost ;         %  Use to calculate logpost
logpostb = logpost ;         %  Use to calculate logpost
trsd = asr ;
ssa = logpost ;                  
sigma2 = logpost ;
	if q==2;
	niter=ones(m,1);
	else;niter=[10;3*ones(m-1,1)];
	end;
id = diag([zeros(1,degree+1)  ones(1,nknots)]) ;

	for i=1:m ;
       	
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
		
		xxid = xx + pwt(i)*id2 ;
		b = inv(xxid) ;
		beta(:,i) = xxid\xy ;
		end ;		
	xxb = xx*b ;
	trsd(i) = trace(xxb) ;
	yhat(:,i) = xm*beta(:,i) ;
	asr(i) = mean( (y-yhat(:,i)).^2 ) ;
	Delta = n + r + 2 + 2*(r1 + nknots + 2)/q ;
	sigma2g = logspace(-5,6,200)' ;
	ssag = n*asr(i) + pwt(i)* norm(beta(degree+2:size(beta,1),i),q).^q  ...
	+ s  + s1 *  (pwt(i)^(2/q))  *  (   sigma2g.^(-2/q+1)   ) ;
	crit =  ssag./(sigma2g) + ...
		Delta*log(sigma2g);
% 	ssag and crit are a function of sigma (column vectors) 
%	crit is the quantity in square brackets in (3) of the paper

	[mincrit,f1] = min(crit') ;  % minimize crit over sigma
	f1=min(f1) ;
%	sigma2, ssa, logpost1 and logspostb are functions of pwt 
	sigma2(i) = sigma2g(f1) ;
	ssa(i) = ssag(f1) ;
	logposta(i)= (-1/2)*crit(f1) ;
%	logposta = 1st 2 lines of (3) in the paper maximized over sigma
		
	
	end ;

eta = 2*log( 2^(1+1/q) *gamma(1/q) / q ) ;

%	logpost is the log posterior as
%	a function of pwt (rows) and s1 (columns)
%
logpost = logposta + (r1+nknots+2)*log(pwt') ...
	 / q - nknots * eta /2 + r1*ones(m,1)*log(s1)/2;

%	logpost is line (3) of the paper

gcv = asr./ ((1-df./n).^2) ;
imax = find(logpost==max(logpost)) ;
imingcv = find(gcv==min(gcv)) ;
yhatgcv = yhat(:,min(imingcv));
yhat = yhat(:,min(imax)) ;
beta = beta(:,min(imax)) ;
sigma2hat = sigma2(min(imax)) ;
df = trsd(min(imax)) ;
