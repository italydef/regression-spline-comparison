function y = pdfmn(x,m,v) ;
% Returns multivariate normal probability density function
% with mean m and covariance matrix v
% Use:  y = pdfmn(x,m,v)
[p1,p2] = size(x);
if p1 == 1; x = x' ; p= p2 ; else ; p = p1 ;end ;
[p1,p2] = size(m) ;
if p1 == 1; m = m' ; end;
y = (2*pi)^(-p/2) .* (det(v)^(-1/2)) .* exp(- .5 * (x-m)' * inv(v).^2 * (x-m) ); 
