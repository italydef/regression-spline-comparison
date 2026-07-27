function fit = LogitReg01(x,y,betahatstart) ;
%  Multivariate logistic regression
%
%  Last edited:  11/5/99
%
%  USAGE:  fit = LogitReg01(x,y,betahatstart) ;
%
%		INPUT (required)
%	x: n by (p-1) matrix of covariates.  Do not include a column of ones;
%		A column of ones is automatically generated.
%	y: n by 1 vector of binary (0-1) responses
%
%		INPUT (optional)
%	betahatstart: starting value of betahat (default is all 0's)
%
%		OUTPUT
%	A structure called 'fit' containing:
%	betahat = MLE of beta
%	yhat = predicted values
%	varmatrix = var/cov matrix of betahat
%	se = standard errors of betahat
%	t = t-values for betahat
%	ci = approximate 95% confidence intervals for betahat
%	betahatseq = sequence of iterates of betahat

if size(y,1) == 1 ;y = y' ; end ;
if size(x,1) < size(x,2) ;x = x' ; end ;
n=rows(x) ;
x = [ones(n,1) x ] ;
p = size(x,2) ;

if nargin < 3 ;
	betahatstart = zeros(p,1) ;
end ;

betahat = betahatstart ;
betahatseq = betahatstart ;
niter = 25 ;
i=1 ;
crit =1 ;
while (i < niter) & (crit > .0000001) ;
	betahatold = betahat ;
	yhat = logistic(x*betahat) ;
	w = yhat .* (1-yhat) ;
	W = kron(ones(1,p),w) ;
	A = x' * (y - yhat) ;
	varmatrix = inv(x'* (x.*W)) ;
	betahat = betahat + varmatrix*A ;
	betahatseq = [betahatseq betahat] ;
	crit = norm(betahat-betahatold,1) ./ (norm(betahatold,1)+eps) ;
end ;
se = sqrt(diag(varmatrix)) ;
t = betahat./se ;
z = 1.96 ;
ci = [betahat - se.*z betahat + se.*z ] ;
fit=struct('betahat',betahat,'yhat',yhat,'varmatrix',varmatrix, ...
	'se',se,'t',t,'ci',ci,'betahatseq',betahatseq) ;

