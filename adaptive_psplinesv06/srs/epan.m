function y = epan(x);
% Computes the epan kernel
% Usage : y = epan(x);
% Returns: y = (3/4).*(x.*x.<1).*(1- x.*x) ) ;
%
i = (x.*x < 1) ;
x = x .*i ;
y = (3/4) * i .* (1- x.*x) ;

