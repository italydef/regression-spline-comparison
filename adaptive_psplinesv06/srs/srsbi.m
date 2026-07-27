function [yhat,beta,mhat,cp,imin,df] = srsbi(x,xgrid,y,string);
%
%   Fits a regression spline to bivariate x's with a quadratic penalty
%
%   USAGE:  [yhat,beta,cp,imin] = srsbi(x,xgrid,y,string)
%
%   Defaults are degree = 2
%                numberknots = max(5,round(.5*sqrt(n)))  
%
n = rows(y) ;
degree = 2 ;                              %  Set defaults
numberknots = max(5,round(.2*sqrt(n))) ;
penaltyweight = logspace(-2,8,25)' ;
eval(string) ;                            %   Change defaults

x1 = x(:,1) ;
x2 = x(:,2) ;
xsort1 = sort(x1) ;                         %   Set knots
loc = n*(1:numberknots)' ./ (numberknots+1) ;
knots1=xsort1(round(loc)) ;

xsort2 = sort(x2) ;                        
loc = n*(1:numberknots)' ./ (numberknots+1) ;
knots2=xsort2(round(loc)) ;
x1=[x1 ; xgrid(:,1)] ;
x2=[x2 ; xgrid(:,2)];


	for i=1:degree ;
		if i == 1 ;
		xm1a = x1 ;
		xm2a = x2 ;
		else ;
		xm1a = [xm1a x1.^i] ;
		xm2a = [xm2a x2.^2] ;
		end ;
	end ;

	for i=1:degree ;
		for j=1:degree;
		xm12a = xm1a(:,i).*xm2a(:,j) ;
		end ;
	end ;
	
	for i=1:(numberknots) ;
		if i == 1 ;
		xm1b = ((x1-knots1(i)).^degree).*(x1 > knots1(i)) ;
		xm2b = ((x2-knots2(i)).^degree).*(x2 > knots2(i)) ;
		else ;
		xm1b = [xm1b ((x1-knots1(i)).^degree).*(x1 > knots1(i))] ;
		xm2b = [xm2b ((x2-knots2(i)).^degree).*(x2 > knots2(i))] ;
		end ;
	end ;

	for i=1:numberknots ;
		for j=1:numberknots ;
			if i==1 & j==1;
			xm12b=xm1b(:,i).*xm2b(:,j) ;
			else ;
			xm12b= [xm12b xm1b(:,i).*xm2b(:,j)] ;
			end ;
		end ;
	end ;

xm = [ones(rows(xm1a),1) xm1a xm2a xm12a xm1b xm2b xm12b ] ;
	

	xx =xm(1:n,:)'*xm(1:n,:) ;
	xy = xm(1:n,:)'*y ;
kk = (1+degree)^2 ;
id = diag([zeros(1,kk) ones(1,cols(xm)-kk)]) ;

m = max(size(penaltyweight)) ;
beta = zeros(cols(xm),m) ;
yhat = zeros(n,m) ;
asr = zeros(m,1) ;
cp = asr ;
trsd = asr ;

	for i=1:m ;
	b = inv(xx + penaltyweight(i)*id) ;
	xxb = xx*b ;
	trsd(i) = trace(xxb) ;
	
	
	beta(:,i) = b * xy ;
	yhat(:,i) = xm(1:n,:)*beta(:,i) ;
	asr(i) = (norm(y-yhat(:,i)).^2)/n ;
		if i==1;
		trsdsd = trace(xxb*xxb) ;
		sigma2hat= n*asr(i)/(n-2*trsd(i)+trsdsd) ;
		end ;
	cp(i) = asr(i) + 2*trsd(i)*sigma2hat/n ;
	
	end ;
imin = find(  (cp==min(cp))  ) 
 

mhat = xm((n+1):(n+rows(xgrid)),:)*beta(:,imin) ;
df = trsd ;
