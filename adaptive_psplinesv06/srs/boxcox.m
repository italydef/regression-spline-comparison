function y = boxcox(x,lambda) ;
%	Usage: y = boxcox(x,lambda)
%
if abs(lambda) < 10*eps;
y=log(x) ;
else;
y=(x.^lambda - 1)./lambda ;
end ;

