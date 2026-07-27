cd('/Volumes/LACIE-WORK 1/Recherche/Anneleen/nnls-anneleen/carroll/adaptive_psplinesv06/srs');

% HeaviSine Wang
n=256;
u=linspace(0, 1, n);
x=u(:);
t=2*sin(4*pi*x)-6*abs(x-0.4).^(3/10)-sign(0.7-x);
SNR=3;
sig=sqrt(var(t));
e = normrnd(0,sig,n,1);
y=t+e;
fit = srslocalauto(x,y);
plot(x,t,'-',x,y,'o',x,fit.yhat(:,1),'b-',x,fit.yhat(:,2),'g-')




setwd('/Volumes/LACIE-WORK 1/Recherche/Anneleen/nnls-anneleen/carroll/adaptive_psplinesv06/srs')
library(R.matlab)
Matlab$startServer()
matlab <- Matlab()
isOpen <- open(matlab)
isOpen
setVariable(matlab, x = xx)
setVariable(matlab, y = as.vector(ydat[1,]))
evaluate(matlab,'fit=srslocalauto(x,y); z=fit.yhat(:,2);')
z<-getVariable(matlab, "z")
yhat=z$z[,1]
