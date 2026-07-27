% 	Program to test latinhc
%
n = 10
p = 1/6
design = latinhc(n,3) ;
plot(design(:,1),design(:,2),'+') ;
phat = mean( (design(:,1) < design(:,2)) & (design(:,2) < design(:,3))  ) 
set(gca,'xlim',[0,1]) ;
set(gca,'ylim',[0,1]) ;

std = sqrt(p*(1-p)/n) 

t = (phat-p)/std
