function y = rounddigit(x,d) ;
% Rounds the elements of x to d digits after the decimal.
% USAGE: y = rounddigit(x,d) ;
%
%	Last edit: 7/21/97
%
y = round(10^d * x)/10^d ;
