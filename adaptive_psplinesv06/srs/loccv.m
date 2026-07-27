function [yhatloccv,spansloccv,hloccv,yhatcv,spancv,hcv] = ...
		loccv(x,y,k,span,kappa) ;
%
% 	Program to compute a local polynomial estimate with a local cv 
%	selected span.  The local span function is a linear spline with
%	"kappa" knots from min(x) to max(x).
%	Also computes a global cv estimate.
%
%		INPUT (Required)
%	x = independent variable (n by 1)
%	y = dependent variable (n by 1)
%
%		INPUT (Optional)
%	k = degree of the local polynomials (default is 1)
%	span = vector of possible values of the global span (this program
%		finds the global span value that minimizes cv) (default of
%		span is linspace((4+2*k)/n,1,20))
%	kappa = number of knots in the local span function (default = 5)
%
%		OUTPUT
%	yhatloccv = estimate at local span that minimizes cv (n by 1)
%	spansloccv = local span that minimizes cv
%	hloccv = bandwidth corresponding to the local span
%	yhatcv = estimate at global span that minimizes cv (n by 1)
%	spancv = global span that minimizes cv
%	hcv = bandwidth corresponding to the local span
%
% 	USAGE: [yhatloccv,spansloccv,hloccv,yhatcv,spancv,hcv] = ...
%		loccv(x,y,k,span,kappa) ;
%
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
if nargin < 6 ;
	kappa = 5 ;
end ;
cv = zeros(length(span),1) ;
for i=1:length(span) ;
	h(:,i) = bspan(x,x,span(i)) ;
	[f(:,i),cv(i)]=lpolycv(x,y,h(:,i),k) ;
end ;
[cvmin,imin]=min(cv) ;
imin = min(imin) ;
spancv = span(imin) ;
hcv= bspan(x,x,spancv) ;
yhatcv= lpolydb(x,y,hcv,x,0,k) ;

knots = linspace(min(x),max(x),kappa)' ;
knots = [min(x)-range(x)/6; knots; max(x)+range(x)/6] ;
heights = spancv*ones(kappa,1) ;
B = bspline(x,kappa-1,1) ;
spangrid = B*heights ;

for iter = 1:5 ;
	for i=2:kappa+1 ;
		induse = (x>knots(i-1)) & (x<knots(i+1)) ;
	
		heightstrial = heights ;
		heightstrial(i-1) = heights(i-1)*(5/4) ;
		spantrial = B*heightstrial ;
		htrial = bspan(x,x,spantrial) ;
		[f,cv1]=lpolycv(x,y,htrial,k,induse) ;
	
		heightstrial = heights ;
		heightstrial(i-1) = heights(i-1)*(4/5) ;
		spantrial = B*heightstrial ;
		htrial = bspan(x,x,spantrial) ;
		[f,cv2]=lpolycv(x,y,htrial,k,induse) ;
	
		heightstrial = heights ;
		spantrial = B*heightstrial ;
		htrial = bspan(x,x,spantrial) ;
		[f,cv3]=lpolycv(x,y,htrial,k,induse) ;
	
		[cvmin,imin] = min([cv1 cv2 cv3]) ;
		if imin == 1 ; 
			heights(i-1) = heights(i-1)*1.2 ;
			elseif imin==2 ;heights(i-1) = heights(i-1)*.8 ;
		end ;
	
	end ;
end ;
spansloccv = B*heights ;
hloccv= bspan(x,x,spansloccv) ;
yhatloccv= lpolydb(x,y,hloccv,x,0,k) ;
