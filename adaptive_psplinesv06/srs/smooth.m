function [Q,B,f]=smooth(x,y,alpha)

n=length(x);

[Q,R,h]=qmtx(x);  % get Q, R, and h as defined in G & S, p12

ydiff=y(2:n)-y(1:n-1);              % cheap way to get Q'*y
Qy=(ydiff(2:n-1)./h(2:n-1)-ydiff(1:n-2)./h(1:n-2));

B=sparse(R+alpha*Q'*Q);

G=chol(B);  %% Note--if chol() is fed a sparse mtx, it
            %% returns a sparse mtx, so in the next 2 steps 
            %% we are still working with a sparse mtx--good.
u=G'\Qy;
v=G\u;

f=y-alpha*Q*v;
