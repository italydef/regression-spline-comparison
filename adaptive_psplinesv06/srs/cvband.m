function [yhat,hmin,cv,spancv] = cvband(x,y,k,span) ;
% Program to compute a local polynomial estimate with a cv selected span
%
%		INPUT (Required)
%	x = independent variable (n by 1)
%	y = dependent variable (n by 1)
%
%		INPUT (Optional)
%	span = vector of possible values of the span (this program
%		finds the value that minimizes cv (default is
%		linspace((4+2*k)/n,1,20))
%	k = degree of the local polynomials (default is 1)
%
%		OUTPUT
%	yhat = estimate at span that minimizes cv (n by 1)
%	hmin = vector of bandwidths corresponding to cv chosen span
%		(n by 1)
%	cv = value of cv as a function of span (dimension = length(span) by 1)
%	spancv = span that minimizes cv
%
% 	USAGE: [yhat,hmin,cv] = cvband(x,y,k,span) ;
%	CALLS: bspan, lpolycv, lpolydb
%
%	Last edit: May 3, 1998
%
%	Copyright: David Ruppert
%
n=length(x) ;
if nargin < 3 ;
k = 1;
end ;
if nargin < 4 ;
span = linspace((4+2*k)/n,1,20) ;
end ;

cv = zeros(length(span),1) ;
for i=1:length(span) ;
h(:,i) = bspan(x,x,span(i)) ;
[f(:,i),cv(i)]=lpolycv(x,y,h(:,i),k) ;
end ;

[cvmin,imin]=min(cv) ;
hmin = h(:,imin) ;
yhat= lpolydb(x,y,hmin,x,0,k) ;
spancv=span(imin) ;
