function [theta0, theta] = MultBC01(X) ;
[n,p]= size(X) ;
if p > n ;
	X = X'
	[n,p] = size(X) ;
end ;
theta0=zeros(p,1) ;
for i=1:p ;
fit = BoxCoxReg01('',X(:,i)) ;
theta0(i) = fit.lambdahat ;
end ;

theta = fminunc('MultBCLik',theta0,[],X) ;
