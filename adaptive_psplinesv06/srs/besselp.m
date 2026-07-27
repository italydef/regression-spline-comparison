%  Computes the bessel function based covariance function in Handcock and 
%  Wallis

n = 50 ;
x = linspace(.001,3,n)' ;
theta1 =  1 ;

y = zeros(n,4) ;
theta2 = [.25; .5; 1; 2] ;
	for i=1:4 ;
	theta1p = theta1 /(2*sqrt(theta2(i))) ;
	y(:,i) = besselk(theta2(i),(x./theta1p)) ;
	y(:,i) = 1 ./ (2^(theta2(i) - 1).*gamma(theta2(i))) ...
	.* ((x./theta1p).^theta2(i)).* y(:,i) ;
	end ;
plot(x,y) ;
