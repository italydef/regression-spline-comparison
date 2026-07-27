function [y,x] = piece(IN,X,Y)
% plots piecewise regression with response Y,
% predictor X, knots IN;

K_all=length(IN);
my_X = [ones(size(X)) ];
for j=1:K_all,
   if IN(j) > 0,
      t=X(IN(j));
      my_X = [my_X max(0,X-t)];
   end;
end;

Beta=my_X\Y;
m=min(X); M=max(X);
x=[m:(M-m)/1000:M]';

my_x = [ones(size(x)) ];
for j=1:K_all,
   if IN(j) > 0,
      t=X(IN(j));
      my_x = [my_x max(0,x-t)];
   end;
end;

y=my_x*Beta;
plot(x,y,'-',X,Y,'o');

