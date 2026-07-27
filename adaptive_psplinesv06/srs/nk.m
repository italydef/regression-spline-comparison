
function Y=nk(p)

Y=zeros(p+1,p+1);
x=linspace(0,p,p+1);
xx=ones(p+1,1)*x;
xx=xx+xx';
Y=muk(xx);


