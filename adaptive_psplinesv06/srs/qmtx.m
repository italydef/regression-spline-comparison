function [SQ,SR,h]=qmtx(x)
% Takes the vector of x values (the explanatory variable)
% and returns the Q mtx as defined
% as defined in Green and Silverman, p12.

n=length(x);

h=x(2:n)-x(1:n-1);
h1=h.^(-1);
h2=-h1(1:n-2)-h1(2:n-1);

Q1=zeros(n,n-2);
Q1(1:n-2,:)=diag(h1(1:n-2));
Q2=zeros(n,n-2);
Q2(3:n,:)=diag(h1(2:n-1));
Q3=zeros(n,n-2);
Q3(2:n-1,:)=diag(h2(1:n-2));
SQ=sparse(Q1)+sparse(Q2)+sparse(Q3);

r1=(h(1:n-2)+h(2:n-1))/3;
r2=h(2:n-2)/6;
SR=sparse(diag(r1))+sparse(diag(r2,-1))+sparse(diag(r2,1));


