function y = huberwt(x,k) ;
%
%
x = x./mad(x) ;
i = (abs(x) < k) ; 
y = i + (1-i).*k./(abs(x)+eps) ;
