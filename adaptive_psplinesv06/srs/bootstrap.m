function y = bootstrap(x);
% Returns a resample of the rows of x;
% Usage:  y = bootstrap(x);
n = rows(x) ;
y = x( ceil(n*rand(n,1)),: ) ;
