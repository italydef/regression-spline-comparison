function y=roundd(x,n) ;
%	Rounds x to n decimal places
%	Example: roundd([.5555 .343],2) = [.56 .34] ;
%	Example: roundd(257,-1) = 260
%	USAGE: y=roundd(x,n) ;
%	Last edit: June 22, 1998
%
y=round(10^n*x)/10^n ;
