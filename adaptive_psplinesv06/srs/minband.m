function minb = minband(xgrid,x,p);
% 
% Computes the minimum bandwidth at each point of xgrid when fitting local
% p degree polynomials with x the vector of x values.
%
%  Last edited:  12/31/95
%
% USAGE: minb = minband(xgrid,x,p);
%
%
n = rows(xgrid);
minb = zeros(n,1) ;
nx = rows(x) ;
x=sort(x) ;

for i=1:nx ;
	if i==1;x2=x(i);
	elseif x(i)>x(i-1); x2 = [x2;x(i)] ;
	end ;
end;

for i = 1:n ;
x3 = sort(abs(x2-xgrid(i))) ;
minb(i) = x3(p+1) ;
end ;
