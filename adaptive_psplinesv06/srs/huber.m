function y = huber(x,k);
% 
% Computes the Huber psi function of x with cutoffs -k and +k
%
i1 = (x<-k) ;
i2 = (x>k) ;
y = -k.*i1 + k.*i2 + x.*(1-i1-i2) ;
