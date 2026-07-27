function y = hubersqr(x,k) ;
% 
y = huber(x,k).^2 -  (  2*(1-cdfn(k))*k*k  +  (2*cdfn(k)-1)  - 2*k*pdfn(k)  ) ;
