function y = huberwt(x,k) ;
%
%
i = (abs(x) < k) 
y = i + (i-1).*k./(abs(x)+eps) ;
