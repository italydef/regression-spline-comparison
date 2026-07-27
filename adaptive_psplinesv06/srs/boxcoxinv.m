function y = boxcoxinv(x,lambda) ;
%	Computes the inverse of the Box-Cox power transformation
%	Usage: y = boxcoxinv(x,lambda)
%
if abs(lambda) < 10*eps;
y=exp(x) ;
else;
y=(1+lambda*x).^(1/lambda);
end ;

