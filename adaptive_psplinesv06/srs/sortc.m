function y = sortc(x,j) ;
%  Sorts the rows of x in ascending order according to the jth column of x
%  Like the GAUSS function sortc
%  USAGE: function y = sortc(x,j) ;
%
[m,n] = size(x) ;
[z,I] = sort(x(:,j)) ;
y=x(I,:) ;

