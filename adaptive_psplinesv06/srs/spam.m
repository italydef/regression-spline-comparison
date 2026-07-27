
function [betahat,mhat,sigma2,var] = spam(x,z,y,h,k) ;

[f,smatrix]=lpolymatrix(z,y,h,k) ;
n = rows(x) ;
proj = (inv(x'*(eye(n)-smatrix)*x))  *  x'*(eye(n)-smatrix) ;
betahat = proj*y ;
var = proj*proj' ;
mhat=smatrix*(y-x*betahat) ;
sigma2 = mean(  (y-x*betahat - mhat).^2  ) ;
