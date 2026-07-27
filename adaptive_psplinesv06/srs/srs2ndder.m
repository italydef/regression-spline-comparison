function [yhat,beta,cp,imin,df] = srs2ndder(x,y,string);
%
%   Fits a regression spline to univariate x's with a quadratic penalty
%
%   USAGE:  [yhat,beta,cp,imin] = srs2ndder(x,y,string)
%
%   Defaults are degree = 2
%                nknots = max(5,round(.5*sqrt(n))) 
%                penwt = logspace(-6,3,30)' 
%
n = rows(x) ;
degree = 2 ;                              %  Set defaults
nknots = max(5,round(.5*sqrt(n))) ;
penwt = logspace(-7,6,30)' ;
w = ones(n,1) ;
eval(string) ;                            %   Change defaults

xsort = sort(x) ;                         %   Set knots
loc = n*(1:nknots)' ./ (nknots+1) ;
knots=xsort(round(loc)) ;
xm=ones(n,1);

	for i=1:degree ;
	xm = [xm x.^i] ;
	end ;
	
	for i=1:(nknots) ;
	xm = [xm ((x-knots(i)).^degree).*(x > knots(i))] ;
	end ;

	xmw = xm.*(w*ones(1,cols(xm))) ;
	yw = y.*w ;
	xx = xmw'*xmw ;
	xy = xmw'*yw ;


m = max(size(penwt)) ;
beta = zeros(cols(xm),m) ;
yhat = zeros(n,m) ;
asr = zeros(m,1) ;
cp = asr ;
trsd = asr ;
knots2 = [min(x); knots; max(x) ] ;

	for i=1:(nknots+1) ;
	 if i==1;
	 M=[2*sqrt(knots2(i+1)-knots2(i))*ones(1,i) zeros(1,(nknots+1-i)) ];
	 else;
         M = [M;2*sqrt(knots2(i+1)-knots2(i))*ones(1,i) zeros(1,(nknots+1-i))] ;
	 end ;
	end ;
	
R = [ones(2,nknots + 3) ;
     ones(nknots+1,2) M'*M ];
     
	for i=1:m ;
	b = inv(xx + penwt(i)*R) ;
	xxb = xx*b ;
	trsd(i) = trace(xxb) ;
	
	
	beta(:,i) = b * xy ;
	yhat(:,i) = xm*beta(:,i) ;
	asr(i) = (norm(y-yhat(:,i)).^2)/n ;
		if i==1;
		trsdsd = trace(xxb*xxb) ;
		sigma2hat= n*asr(i)/(n-2*trsd(i)+trsdsd) ;
		end ;
	cp(i) = asr(i) + 2*trsd(i)*sigma2hat/n ;
	
	end ;
imin = find(  (cp==min(cp))  ) ;
df= trsd ;
