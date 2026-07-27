function [x2,y2,varyfactor,varest] = prebin2(x,y,nbin,deg,order) ;
%
%      copied from prebin
%
%      Last edited 6/16/95
%
%  Bins (x,y) data using local polynomial regression with "boxcar" weights
%
%  The (x,y) data are sorted according to their x-values and put into
%  nonoverlapping groups of approximately "nx" data points.  For the ith
%  group 
%    x2(i) = mean of the x's
%    y2(i) = estimated derivative of the regression function
%    varest(i) = MSE from this fit
%    varyfactor = factor to multiply by varest(i) to get the estimated
%                 variance of y2(i)
%   
%
%   INPUT: x,y = raw data
%          nbin = number of bins
%          deg = degree of the local polynomials
%          order = order of the derivative returned
%
%
%   The program computes n=dim(x) and nx= floor(n/nbin)
%
%   USAGE: function [x2,y2,varyfactor,varest] = prebin(x,y,nbin,deg,order) ;
% 
%
global nx ;
[x,I]=sort(x) ;
y=y(I) ;
n= rows(x) ;
nx = floor(n/nbin) ; 

x2=zeros(nbin,1) ;
y2 = x2 ;
varyfactor = x2 ;
varest = x2 ;

	for m = 0:(nbin-1) ;
            if m < nbin-1 ;
	    x3 = [ ones(nx,1) x((m*nx+1):((m+1)*nx)) ] ;
	    y3 = y((m*nx+1):((m+1)*nx)) ;
	    else;
	    x3=[ ones(n-m*nx,1) x((m*nx+1:n)) ] ;
	    y3 = y((m*nx+1):n) ;
	    end ;
	xmean = mean(x3(:,2)) ;
	x2(m+1) = xmean ;
	x3(:,2) = x3(:,2) - xmean ;
		for i = 2:deg ;
		x3 = [x3 x3(:,2).^i] ;
		end ;
	ixx = inv(x3'*x3) ;
	beta = ixx*(x3'*y3) ;
	y2(m+1) = fact(order)*beta(1+order) ;
	varest(m+1) =sum((y3-x3*beta).^2)/(rows(x3)-(1+deg)) ;
	varyfactor(m+1) = (fact(order).^2).*ixx(1+order,1+order) ;
	end ;


