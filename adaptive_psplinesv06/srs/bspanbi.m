function bsp = bspanbi(xgrid,x,p);
% 
%   Computes the bandwidth at each point of xgrid corresponding to
%    a span of p with x the vector of x values for bivariate x and xgrid
%
%  Last edited:  7/5/96
%
% USAGE: bsp = bspanbi(xgrid,x,p);
%
%
n = rows(xgrid);
nx = rows(x) ;
bsp = zeros(n,1) ;
np = ceil(nx*p) ;
for i = 1:n ;
x2 = sort(     sqrt(    (  x-ones(nx,1)*xgrid(i,:)  ).^2 *ones(2,1)    )     ) ;
bsp(i) = x2(np) ;
end ;
