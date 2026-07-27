function fit=BoxCoxReg01(x,y,param) ;

if nargin < 3 ;
	param = '' ;
end ;

if isfield(param,'lamL') == 0 ;
	lamL = -1.5 ;
else ;
	lamU = param.lamL ;
end ;

if isfield(param,'lamU') == 0 ;
	lamU = 1.5 ;
else ;
	lamU = param.lamU ;
end ;

if isfield(param,'Nlambda') == 0 ;
	Nlambda = 200 ;
else ;
	Nlambda = param.Nlambda ;
end ;
n=length(y) ;

if isempty(x) == 1 ;
	x = ones(n,1) ;
else ;
	x = [ones(n,1) x] ;
end ;

lambda = linspace(lamL,lamU,Nlambda)' ;
loglik = 0*lambda ;
logy = sum(log(y)) ;
xx = x'*x ;

for i=1:Nlambda ;
ylam = boxcox(y,lambda(i)) ;
beta = xx\(x'*ylam) ;
sig2 = mean((ylam-x*beta).^2) ;
loglik(i) = -.5*n*log(sig2) + (lambda(i)-1)*logy ;
end ;	%	End "for i=1:Nlambda" 

imax = min(find(loglik == max(loglik))) ;
lambdahat = lambda(imax) ;
ylam = boxcox(y,lambdahat) ;
betahat = xx\(x'*ylam) ;
sigma2hat = mean((ylam-x*betahat).^2) ;
fit = struct('lambdahat',lambdahat,'betahat',betahat,'sigma2hat',sigma2hat) ;
