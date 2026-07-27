function grid = addgrid(x,j,string) ;
eval(string) ;
n = rows(x) ;
m = cols(x) ;
	for i=1:m;
		if i==j;
		xcol=linspace(min(x(:,j)),max(x(:,j)),ngrd) ;
		else ;
		xcol = ones(ngrd,1)*mean(x(:,j)) ;
		end ;
	
		if i==1 ;
		xnew = xcol ;
		else;
		xnew = [xnew xcol] ;
		end ;
	end ;	
xold = x ;
x = xnew ;

d = cols(x) ;
xm1=ones(ngrd^2,1);
	for j = 1:d ;
	xj = x(:,j) ;                 %  CREATE KNOTS AT
	xsort = sort(xold(:,j)) ;            %  SAMPLE QUANTILES
	loc = n*(1:nknots+2)' ./ (nknots+3) ;
	knots=xsort(round(loc)) ;
	knots=knots(2:nknots+1) ;  %  REMOVE KNOTS NEAR BOUNDARIRES FOR
                                        %  STABILITY (= LOW VARIABILITY)
		for i=1:degree ;       
			if i==1 ;xmj = xj ;
			else ;
			xmj = [xmj xj.^i] ;
			end ;
		end ;
		
		for i=1:nknots ;
		xmj = [xmj ((xj-knots(i)).^degree).*(xj > knots(i))] ;
		end ;
	xm1 = [xm1 xmj] ; 
	end ;

xm2=ones(ngrd^2,1);
	for j = 1:d ;
	xj = x(:,j) ;                 %  CREATE KNOTS AT
	xsort = sort(xold(:,j)) ;            %  SAMPLE QUANTILES
	loc = n*(1:nknots2+2)' ./ (nknots2+3) ;
	knots2=xsort(round(loc)) ;
	knots2=knots2(2:nknots2+1) ;  %  REMOVE KNOTS NEAR BOUNDARIRES FOR
                                        %  STABILITY (= LOW VARIABILITY)
		for i=1:degree ;       
			if i==1 ;xmj = xj ;
			else ;
			xmj = [xmj xj.^i] ;
			end ;
		end ;
		
		for i=1:nknots2 ;
		xmj = [xmj ((xj-knots(i)).^degree).*(xj > knots(i))] ;
		end ;
	xm2 = [xm2 xmj] ;
	end ;
xm2=xm2(:,2:cols(xm2)) ;

degnk2 = degree + nknots2 ;
xm3 = ones(ngrd^2,1) ;
size(xm2)
degnk2
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
grid = xm ;
