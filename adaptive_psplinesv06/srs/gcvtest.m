function gcv=gcvtest(Q,B,ss,n,alpha)

%% This returns the GCV score for a given alpha 
%% where f is the vector of fitted values using 
%% cubic smoothing splines, Q & B (B=R+alpah*Q'*Q) as defined
%% in Green and Silverman (p.34 bottom).  In this function, let A
%% represent the matrix A(alpha)=I-alpha*Q*inv(B)*Q'
%% as defined on p.34. 


Binv=inv(B);
A=eye(n)-alpha*Q*Binv*Q';

den=(1-trace(A)/n)^2;
gcv=(ss/den)/n;


