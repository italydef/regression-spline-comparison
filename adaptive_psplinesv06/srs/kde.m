function y = kde(x,grid,h) ;
% Return a kernel density estimate of x, on grid, with bandwidth std(x)*h.
% Usage:  y = kde(x,grid,h) ;
n = size(x,1) ;
m = size(grid,1) ;
y = zeros(m,1) ;
h = h * std(x) ;
for i = 1:m ;
y(i) = sum( epan((x-grid(i))./ h) ) ./ (n*h) ;
end ;
