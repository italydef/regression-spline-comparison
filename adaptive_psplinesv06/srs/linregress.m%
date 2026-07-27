function   fit = linregress(x,y,intercept,w) ;
%
%	Subroutine for linear regression analysis
%
%		INPUT (required)
%	x= matrix of predictor variables
%	y=vector of responses
%
%		INPUT (optional)
%	intercept (1 => use an intercept (default), ~1 => no intercept
%	w = weights (for weighted least squares so that 
%	\sum (y_i - x_i\beta)^2 w_i is minimized
%		is 1)
%
%		OUTPUT (In a structure called "fit") 
%	beta = least squares estimate
%	se = standard errors of beta
%	mse = mean square error
%	yhhat = predicted values
%	ssr = regression sum of squares
%	sse = error sum of squares
%	covbeta = estimated covariance matrix of beta
%	resid = raw residuals
%	stdres = standardized residuals
%	hi = hat diagonals
%	press = Press statistics
%	t = t-statistics for beta = beta ./ se
%	
%
%	USAGE: fit = linregress(x,y,intercept,w) ;
%
%	Last edit:  Oct 8, 1998
%
[n,n2] = size(x) ;

if nargin < 4 ;
	w = ones(n,1) ;
end ;

if nargin < 3 ;
	intercept = 1 ;
end ;

n3 = n2 + 1 ;           % number of parameters including intercept
if intercept == 1 ;
	x = [ones(n,1) x] ;    % add column of ones for the intercept
end ;
xw = x.*(w*ones(1,size(x,2))) ;

xx = xw'*x ;
xy = xw'*y ;
beta = xx\xy ;         % solves the equation xx*beta = xy
yhat = x*beta ;
resid = y - yhat ;
sse = sum( resid.^2 ) ;
ssr = sum( (yhat-mean(y)).^2 ) ;
r2 = ssr/(sse+ssr) ;
mse = sse/(n-n3) ;
covbeta = mse*inv(xx) ;

se = sqrt(diag(covbeta)) ;  % "diag" extracts the diagonal of a matrix
t = beta ./ se ;
p = 2*(1-tcdf(abs(t),n-n3)) ;

hi = diag(x*inv(xx)*x') ;
stdres = resid./sqrt( mse*(ones(n,1)- hi) ) ;
presid = resid./(1-hi) ;
msei = mse .* ( (n-n3-stdres.^2)./ (n-n3-1) ) ;

stdpresid = presid ./ sqrt( msei ./ (1-hi) )  ;

press = sum( presid.^2 ) ;
gcv = sse ./ (n * (1-length(beta)/n)^2 ) ;

rstudent = stdres  .*  sqrt(  (n-n3-1) ./ (n - n3 - stdres.^2) ) ;

fit = struct('beta',beta,'se',se,'mse',mse,'yhat',yhat,'ssr',ssr,...
	'sse',sse,'r2',r2,'covbeta',covbeta,'resid',resid, ...
	'stdres',stdres,'hi',hi,'press',press,'t',t,'p',p,'gcv',gcv, ...
	'rstudent',rstudent,'presid',presid,'stdpresid',stdpresid, ...
	'msei',msei) ;
