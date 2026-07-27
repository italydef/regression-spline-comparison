function [IN,OUT,K] = build_OUT(X)
% constructs the set of potential knots;

[N,p]=size(X);
alpha=.01;
skip=ceil((-log2(-(1/N)*log(1-alpha)))/2.5);
k=1;
OUT(k)=1; 
sofar=0;
for i=2:N,
   if (X(i)~=X(i-1))
      sofar=sofar+1;   
         if sofar==skip
            k=k+1;
            OUT(k)=i;
            sofar=0; 
         end;
   end;
end;
for i=1:k,
   IN(i)=0;
end;
K=k;
