function [fmin,x_sorted,fmin_sorted]=spline_hogan(x,y,alpha)

% Pre:
% Takes the vector of y values (the observed responses), and 
% the vector of x values (the explanatory variable)
% and returns the vector fmin of fitted values using fixed 
% knots natural cubic smoothing splines and GCV to choose 
% smoothing parameter.  Alpha is the grid of values for the
% grid search.
%
%		INPUT (required)
%	x = independent variable
%	y = dependent variable
%
%		INPUT (optional)
%	alpha = grid of values for the grid search (default is
%		alpha = logspace(-8,8,50) ;
%
%		OUTPUT:
%	fmin = fitted values
%	x_sorted = x sorted from smallest to largest
%	fmin_sorted = fmin sorted according to x
%
%	x_sorted and fmin_sorted are useful for plotting
%

% Note: the vector A in this program allows you to track the 
% GCV scores corresponding to the grid of alpha values.  The first 
% has the alpha value used in smoothing, and the associated GCV 
% score is in the second column. 
%
%	Written by Steve Hogan, with minor additions by David Ruppert
%
%	Last edit:	May 28, 1998
%
%	Calls: smooth, gcvtest 
%

if nargin < 3 ;
alpha = logspace(-8,8,50) ;
end ;


n=length(x);
m=length(alpha);
gcv=1e10;

data = [x y (1:n)'] ;
data = sortrows(data,1) ;
x=data(:,1) ;
y=data(:,2) ;

fmin=[];

A=zeros(m,2);

for i=1:m

[Q,B,f]=smooth(x,y,alpha(i));

diff=(f-y);
ss=diff'*diff;

A(i,1:2)=[alpha(i),gcvtest(Q,B,ss,n,alpha(i))];

 if A(i,2)<gcv
    fmin=f;
    gcv=A(i,2);
    aopt=A(i,1);
 end;

end;
x_sorted = x ;
fmin_sorted = fmin ;
data = [data fmin] ;
data = sortrows(data,3) ;
fmin = data(:,4) ;






