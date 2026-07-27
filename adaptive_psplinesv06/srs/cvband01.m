function fit = cvband(x,y,pmean,span) ;
%
% Program to compute a local polynomial estimate with a cv selected span
%
%		INPUT (Required)
%	x = independent variable (n by 1)
%	y = dependent variable (n by 1)
%
%		INPUT (Optional)
%	span = vector of possible values of the span (this program
%		finds the value that minimizes cv (default is
%		linspace((4+2*pmean)/n,1,20))
%	pmean = degree of the local polynomials (default is 1)
%
%		OUTPUT (in a structure called "fit")
%	yhat = estimate at span that minimizes cv (n by 1)
%	hmin = vector of bandwidths corresponding to cv chosen span
%		(n by 1)
%	cv = value of cv as a function of span (dimension = length(span) by 1)
%	spancv = span that minimizes cv
%	xgrid = grid of 100 points from min(x) to max(x)
%	mhat = estimated mean function on xgrid (same estimate as
%		yhat but calculated on xgrid, not the values of x)
%
% 	USAGE: [yhat,hmin,cv] = cvband(x,y,pmean,span) ;
%	CALLS: bspan, lpolycv, lpolydb
%
%	Last edit: June 13, 2000
%
%	Copyright: David Ruppert
%
n=length(x) ;

if nargin < 3 ;
	pmean = 1;
end ;

if nargin < 4 ;
	span = linspace((4+2*pmean)/n,1,20) ;
end ;

cv = zeros(length(span),1) ;

for i=1:length(span) ;
	h(:,i) = bspan(x,x,span(i)) ;
	[f(:,i),cv(i)] = lpolycv(x,y,h(:,i),pmean) ;
end ;

[cvmin,imin]=min(cv) ;
hmin = h(:,imin) ;
yhat = lpolydb(x,y,hmin,x,0,pmean) ;
spancv=span(imin) ;
minx = min(x) ;
maxx = max(x) ;
xgrid = linspace(minx,maxx,100)' ;
band = interp1(x,hmin,xgrid,'linear') ;
mhat = lpolydb(x,y,band,xgrid,0,pmean) ;


fit = struct('yhat',yhat,'hmin',hmin,'cv',cv,'spancv',spancv, ...
	'xgrid',xgrid,'mhat',mhat) ;
