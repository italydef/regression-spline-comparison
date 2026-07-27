% Returns the t-distribution probability density function of x with k 
% degrees of freedom.
% USE:  y = pdft(x,k)
%
function y = pdft(x,k);
y= gamma((k+1)/2) / (sqrt(pi*k) * gamma(k/2));
y= y .* (1+(x.^2)/k).^(-(k+1)/2);
