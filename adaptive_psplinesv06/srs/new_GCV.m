function fun_val = new_GCV(i,IN,OUT,X,Y)
% builds the model with the knots currently in IN
% minus the trial knot, the i-th. Returns the GCV.

if IN(i)==0 
   fun_val=inf;
   return;
end;
K_all=length(OUT);
K=0;
my_X = [ones(size(X)) ];
for l=1:K_all,
   if (IN(l) > 0)&(l~=i),      
      K=K+1;
      t=X(IN(l));
      my_X = [my_X max(0,X-t)];
   end;
end;

[N,p]=size(X);
Beta=my_X\Y;
res=Y-my_X*Beta;
GCV=res'*res/(1-(3*K+1)/N).^2;
fun_val=GCV;


