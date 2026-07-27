function y = mad(x);
%
%   Computes the median absolute deviation from the median
%   USAGE: y = mad(x)
%
%
y = median(abs(x-median(x))) ;
