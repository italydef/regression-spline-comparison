function y = geomean(x);
%  Returns the geometric mean of x
%  USE:   y = geomean(x)
y = exp(mean(log(x)));
