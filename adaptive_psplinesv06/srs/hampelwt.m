function y = hampelwt(x,k1,k2,k3) ;
%
%
const = .674 ;
madx = mad(x)/const ;
x = abs(x./madx) ;
i1 = (x < k1) ; 
i2 = (x < k2) ;
i3 = (x < k3) ;
y = i1 + (1-i1).*i2.*k1./(x+eps) +  i3.*(1-i2).*(k3-x)./(x.*(k3-k2)) ;
