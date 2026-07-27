function gcv=gcvtest(Q,B,ss,n,alpha)

%% This returns the GCV score for a given alpha 
%% where f is the vector of fitted values using 
%% cubic smoothing splines, Q & B (p.34 bottom) as defined
%% in Green and Silverman.  In this function, let A
%% represent the matrix A(alpha)=I-alpha*Q*inv(R+a*Q'*Q)*Q';
%% as defined on p.34.  Let the vector a be the diagonals of 
%% I-A(alpha).  

Binv=inv(B);

a=zeros(n,1);

a(1) = Q(1,1)^2*Binv(1,1);
a(2) = Q(2,1)^2*Binv(1,1)+Q(2,2)^2*Binv(2,2)+2*Q(2,1)*Q(2,2)*Binv(1,2);
a(n-1) = Q(n,n-3)^2*Binv(n-3,n-3) + Q(n-1,n-2)^2*Binv(n-3,n-3) + 2*Q(n-1,n-3)*Q(n-1,n-2)*Binv(n-3,n-2);
a(n) = Q(n,n-2)^2*Binv(n-2,n-2);

for i=3:n-2
   a(i)=Q(i,i-2)^2*Binv(i-2,i-2)+Q(i,i-1)^2*Binv(i-1,i-1)+Q(i,i)^2*Binv(i,i) + 2*Q(i,i-2)*Q(i,i-1)*Binv(i-2,i-1)+2*Q(i,i-2)*Q(i,i)*Binv(i-2,i) + 2*Q(i,i-1)*Q(i,i)*Binv(i-1,i);
end;   

trc=n-alpha*sum(a);
den=(1-trc/n)^2;
gcv=(ss/den)/n;





