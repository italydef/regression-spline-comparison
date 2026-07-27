function y = cdfn(x) ;
% 
% Computes the standard normal cdf using MATLAB erf and erfc functions
%
%  USAGE:  y = cdfn(x) ;
%
%
global i1 i2 ;
i1 = (x > 0 ) ;
i2 = 1 - i1 ;
x = x ./ sqrt(2) ;
y = i1.*(  erf(x) +  1 ) ./ 2+ ...
    i2.* (  erfc(-x)./2  ) ;
