function y=constants(p)

% From Jean Opsomer on 9/5/96
% Function generates the 6 kernel-dependent constants for the
% plug-in bandwidths:
%                 C1(K,p) for h,
%                 C21(K,p), C22(K,p) for g
%                 C3(K,p) for kappa
% for each element of p
%
% Syntax: y=constants(p)
%
% Uses: nk, muk, rk

[pmax,ii]=max(p);
D=max(size(p));
y=zeros(6,D);
yy=zeros(6,1);
for i=1:D
Nmat=inv(nk(p(i)));
x=linspace(p(i)+1,2*p(i)+1,p(i)+1);
mu1=Nmat(1,:)*muk(x)';
x=linspace(0,p(i),p(i)+1);
xx=ones(p(i)+1,1)*x;
xx=xx+xx';
xx=rk(xx);
R1=sum(sum(Nmat(1,:)'*Nmat(1,:).*xx));
yy(1)=((p(i)+1)*fact(p(i))^2*R1/2/mu1^2)^(1/(2*p(i)+3));

Nmat=inv(nk(p(i)+2));
x=linspace(p(i)+3,2*p(i)+5,p(i)+3);
f=fact(p(i)+1);
mu2=f*Nmat(p(i)+2,:)*muk(x)';
Nmat=inv(nk(p(i)+2));
x=linspace(0,p(i)+2,p(i)+3);
xx=ones(p(i)+3,1)*x;
xx=xx+xx';
xx=rk(xx);
R2=f^2*sum(sum(Nmat(p(i)+2,:)'*Nmat(p(i)+2,:).*xx));
temp=(fact(p(i)+3)*R2/2/mu2);
yy(2)=temp^(1/(2*p(i)+5));
yy(3)=(temp*(2*p(i)+3)/2)^(1/(2*p(i)+5));

temp=R1-2*Nmat(1,1)*0.75;
temp1=(1+(temp>0)*(1/(2*p(i)+2)-1))*abs(temp);
yy(4)=(temp1*f^2/mu1^2)^(1/(2*p(i)+3));

yy(5)=(abs(temp)*f^2/mu1^2)^(1/(2*p(i)+3));
yy(6)=temp;
y(:,i)=yy;

end


