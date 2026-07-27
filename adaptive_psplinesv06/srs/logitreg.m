function [beta,yhat,varmatrix] = logitreg(x,y,betastart,weights) ;
%  Univariate logistic polynomial regression
%
%  Last edited:  9/7/98
%
%  USAGE:  function [beta,yhat] = logitreg(x,y,betastart,weights) ;
%
%


p=max(size(betastart)) - 1 ;
if size(y,1) == 1 ;y = y' ; end ;
if size(x,1) == 1 ;x = x' ; end ;
if size(weights,1) == 1 ;weights = weights' ;end ;

n=rows(x) ;
X = ones(n,1) ;
	for i=1:p ;
	X = [X (x.^i)] ;
	end ;

beta = betastart ;
niter = 15 ;
i=1 ;
crit =1 ;
	while (i < niter) & (crit > .000001) ;
	betaold = beta ;
	yhat = logistic(X*beta) ;
	Lprime = yhat .* (1-yhat) ;
	L2 = kron( ones(1,p+1),(weights.*Lprime) ) ;
	l1 = X' * ( weights.*(y - yhat) ) ;
	D = X'* (X.*L2) ;
	Dinv = inv(D) ;
	beta = beta + Dinv*l1 ;
	crit = norm(beta-betaold,1) ./ (norm(betaold,1)+eps) ;
	end ;
L3 = kron( ones(1,p+1),(Lprime.*weights.*weights) ) ;
C = X'*(L3.*X) ;

varmatrix = Dinv*C*Dinv;

