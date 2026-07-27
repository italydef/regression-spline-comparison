% The macro kernel.m computes the the two parts of the product 
%kernel called Kh1 and Kh2. Remember that we only have to compute the
%value of Kh1 and Kh2 in 2(L1+L2+2) points in order to derive the 
%value of the kernel in all gridpoints. The values of Kh1 and Kh2 are 
%then stored in two vectors called Kh1 and Kh2.

function [Kh1,Kh2,L1,L2]=kernel(tau,h1,h2,M1,M2,de1,de2)

  L1=round(min(tau*h1/de1,M1));
  L2=round(min(tau*h2/de2,M2));
         
  for l1=-L1:1:L1
  
   Kh1(l1+L1+1)=exp(-(l1*de1)^2/(2*h1^2));
    
  end

  for l2=-L2:1:L2
  
   Kh2(l2+L2+1)=exp(-(l2*de2)^2/(2*h2^2));
    
  end


