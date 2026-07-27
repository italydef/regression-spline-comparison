function fun_val = new_ASR(i,IN,OUT,X,Y)
% builds the model with the knots currently in IN
% plus the trial knot, the i-th. Returns the ASR.

K_all=length(OUT);
my_X = [ones(size(X)) ];
for j=1:K_all,
   if IN(j) > 0,
      t=X(IN(j));
      my_X = [my_X max(0,X-t)];
   end;
end;
t=X(OUT(i));
my_X = [my_X max(0,X-t)];

Beta=my_X\Y;
res=Y-my_X*Beta;
ASR=res'*res;
fun_val=ASR;


