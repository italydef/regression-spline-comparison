function [yhat,beta,mhat,gcv,imin,pw,df] = ...
	srsinteract(x,y,degree,nknots1,nknots2,penwt1,penwt2);
%           Smoothed Regression Splines for Interaction Model with 
%           with a quadratic penalty.  %
%
%   	Last edited: 12/1/97
%
%	Calls: powerbasis, quantileknots
%
%
if nargin < 3 ;
degree = 2 ;   
end ;

if nargin < 4 ;
nknots1 = 5 ;
end ;

if nargin < 5 ;
nknots2 = 2 ;
end ;

if nargin < 6 ;
penwt1 = logspace(-10,10,20)' ;
end ;

if nargin < 7 ;
penwt2 = logspace(0,20,20)' ;
end;

npw1 = length(penwt1) ;
npw2 = length(penwt2) ;

pw = [kron(ones(npw2,1),penwt1) kron(penwt2,ones(npw1,1))] ;

[n,d] = size(x) ;

xm1 = ones(n,1) ;
xm2 = [] ;
	for j = 1:d ;
	xj = x(:,j) ;                 %  CREATE KNOTS AT
	knots1 = quantileknots(x,nknots1) ;
	knots2 = quantileknots(x,nknots2) ;
	xmj1 = powerbasis(xj,degree,knots1) ;
	xmj2 = powerbasis(xj,degree,knots2) ;
	xm1 = [xm1 xmj1(:,2:size(xmj1,2))] ;
	xm2 = [xm2 xmj2(:,2:size(xmj2,2))] ;

	end ;

degnk2 = degree + nknots2 ;
xm3 = ones(n,1) ;

	for i=1:d-1;
		for j=i+1:d ;
			for k = 1:(degnk2) ;
			for l = 1:(degnk2) ;
			xm3=[xm3 xm2(:,(i-1)*degnk2+k).*xm2(:,(j-1).*degnk2+l)];
			end ;
			end ;
		end ;
	end ;

xm3 = xm3(:,2:cols(xm3)) ;
xm = [xm1 xm3] ;
xx =xm'*xm ;
xy = xm'*y ;

m = max(size(pw)) ;
beta = zeros(cols(xm),m) ;
yhat = zeros(n,m) ;
asr = zeros(m,1) ;
gcv = asr ;
trsd = asr ;


	for i=1:m ;        %  Compute the regression spline for the
                           %  various penalty weights.

	patt = kron(ones(1,d),[zeros(1,degree)  ones(1,nknots1)]) ;
        id = diag([0 pw(i,1)*patt pw(i,2)*ones(1,cols(xm3))]) ;

	b = inv(xx + id) ;
	beta(:,i) = b*xy ;
	xxb = xx*b ;
	trsd(i) = trace(xxb) ;

	yhat(:,i) = xm*beta(:,i) ;
	asr(i) = (norm(y-yhat(:,i)).^2)/n ;   %  asr = average squared residual
		if i==1;
		trsdsd = trace(xxb*xxb) ;
		sigma2hat= n*asr(i)/(n-2*trsd(i)+trsdsd) ;
		end ;
	gcv(i) = asr(i) / (1-trsd(i)/n)^2 ;
	end ;

imin = find(  (gcv==min(gcv))  ) ;
imin = min(imin) ;
yhat = yhat(:,imin) ;
beta = beta(:,imin) ;

for j=1:d ;
xmj = xm(:, 2+(j-1)*(degree+nknots1) : 1+j*(degree+nknots1) ) ;
mhat(:,j)=xmj*beta(2+(j-1)*(degree+nknots1) : 1+j*(degree+nknots1)) ;
end ;
mhat = mhat - ones(n,1)*mean(mhat) ;
df=trsd ;
