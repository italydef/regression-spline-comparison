function y = tent(x,a,b) ;
% 
%  Computes the piecewise linear function whose graph is an isosceles triangle
%  with base the interval (a,b) and whose height is 1.
%
%global a b m c i1 i2 i3 ;
m = (a+b)/2 ; 
c = (b-a)/2 ;
i1 = (x>a) ;
i2 = (x>m) ;
i3 = (x>b) ;
y= (x-a).*i1.*(1-i2)/c   +    (b-x).*i2.*(1-i3)/c ;


