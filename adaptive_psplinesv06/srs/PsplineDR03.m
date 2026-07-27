function  fit = PsplineDR03(x,y,param);
%
%	Fits a P-spline to univariate x's with Demmler-Reinsch
%	algorithm.
%
%	Copied from PsolineDR02 - like that program but optional input
%		enters as a structure.
%
%	Also allows a quadratic integral penalty on the 2nd derivative
%
%
%	USAGE:  fit = PsplineDR03(x,y,param);
%
%
%
%		INPUT - REQUIRED
%	x = independent variable (univariate)
%	y = response (same length as x)
%
%		INPUT - OPTIONAL (put in a structure that is the third
%			argument in)
%	degree = degree of the spline (default is 2) (changed to 3 if
%			smooth_spline_penalty = 1)
%	nknots = number of knots (default min of floor(.3*n) and 20)
%	extrapen: if 1 then the x^degree term is penalized; if 0 then not
%			(default is 0)
%	penwt = trial values of the penalty weight (one is chosen by 
%		minimizing gcv (default is logspace(-12,12,100))
%	boundstab = parameter passed to quantileknots (see that program)
%	knots (default is to generate nknots using the program
%		quantileknots)
%	istd = 1 if x and y are to be standardized before computation
%		(default = 1)
%	smooth_spline_penalty = if 1 then penalty is on the
%		integral from min(x) to max(x) of the square
%		of the second derivative
%	
%		OUTPUT
%	Returns a structure "fit" with the following components:
%
%	yhat = fitted values (same length as x)
%	beta = spline coefficients
%	gcv = vector of gcv values (same length as penwt)
%	imin = index of penwt where minimum of gcv occurs
%	dffit = degrees of freedom of the smoother (same length as penwt)
%	knots = knots (length = nknots)
%	postvarbeta = posterior variance of beta (from a Bayesian model)
%	postvaryhat = posterior variance of yhat (from a Bayesian model)
%	xm = design matrix of the spline (not returned if n > 1000)
%	xx = xm' * xm
%	a =  penwt(imin)
%	penwt
%	sigma2hat
%	dfres
%	varbeta
%	stdyhat
%	yhatder
%	postvaryhatder
%	ulimit
%	llimit
%	ulimitder
%	llimitder
%
%		That is, 
%	fit = struct( ...
%	'yhat',yhat,'beta',beta,'gcv',gcv,'imin',imin,'dffit',dffit, ...
%	'knots',knots,'postvarbeta',postvarbeta,'postvaryhat',postvaryhat, ...
%	'xm',xm,'xx',xx,'a',a,'penwt',penwt,'sigma2hat',sigma2hat, ...
%	'dfres',dfres,'varbeta',varbeta,'stdyhat',stdyhat, ...
%	'yhatder',yhatder,'postvaryhatder',postvaryhatder, ...
%	'ulimit',ulimit,'llimit',llimit,'ulimitder',ulimitder, ...
%	'llimitder',llimitder,'asr',asr) ;
%
%
%	CALLS: powerbasis01, quantileknots
%
%	Last edit: 8/21/2001
%


n = size(x,1) ;

if nargin < 3 ;
	param = '' ;
end ;



if isfield(param,'gcvfact') == 0 ;
	gcvfact = 1 ;
else ;
	gcvfact = param.gcvfact ;
end ;



if isfield(param,'istd') == 0 ;
	istd = 1 ;
else ;
	istd = param.istd ;
end ;


if isfield(param,'boundstab') == 0 ;
	boundstab = 0 ;
else ;
	boundstab = param.boundstab ;
end ;

if isfield(param,'penwt') == 0 ;
	penwt = logspace(-12,12,100) ;
else ;
	penwt = param.penwt ;
end ;

if isfield(param,'extrapen') == 0 ;
	extrapen = 0 ;
else ;
	extrapen = param.extrapen ;
end ;

if isfield(param,'nknots') == 0 ;
	nknots = min([floor(.3*n), 20]) ;
else ;
	nknots = param.nknots ;
end ;


if isfield(param,'knots') == 0 ;
	knots = quantileknots(x,nknots,boundstab) ;
else ;
	knots = param.knots ;
	nknots = length(knots) ;
end ;

if isfield(param,'degree') == 0 ;
	degree = 2 ;
else ;
	degree = param.degree ;
end ;

if isfield(param,'smooth_spline_penalty') == 0 ;
	smooth_spline_penalty = 0 ;
else ;
	smooth_spline_penalty = param.smooth_spline_penalty ;
end ;
	



if istd == 1 ;
	stdx = std(x) ;
	meanx = mean(x) ;
	x = (x - meanx) ./ stdx ;
	meany = mean(y) ;
	y = y - meany ;
	knots = (knots - meanx) ./ stdx ;
end ;


n = length(x) ;

if smooth_spline_penalty == 1 ;
	degree = 3 ;
end ;	%	end "if smooth_spline_penalty == 1 ;"


xm = powerbasis01(x,degree,knots) ;
xx = xm'*xm ;

if extrapen == 1 ;
	id = [zeros(1,degree)  ones(1,nknots+1)] ;
else ;
	id = [zeros(1,degree+1)  ones(1,nknots)] ;
end ;

D = diag(id) ;

if smooth_spline_penalty == 1 ;
maxx = max(x) ;
knots2 = [min(x);knots] ;

for i = 1:4+nknots ;

if i < 3 ; 
D(i,:) = zeros(1,4+nknots) ;
elseif i== 3 ;
D(i,1) = 0 ;
D(i,2) = 0 ;
D(i,3) = 4 ;

for i2 = 4:4+nknots ; 
D(i,i2) = 6*(maxx-knots2(i2-3))^2 ; 
end ;	%	end "for i2 = 4:4+nknots ; "

elseif i > 3 ;
D(i,1) = 0 ;
D(i,2) = 0 ;
D(i,3) = 6*(maxx-knots2(i-3))^2 ;

for i2 = 4:4+ nknots ; 
knotmax = knots2( max([i-3 i2-3]) ) ;
knotmin = knots2( min([i-3 i2-3]) ) ;
D(i,i2) = 36*( (maxx-knotmax)^3/3 + ...
	(knotmax - knotmin) * (maxx-knotmax)^2/2) ;
end ;	%	end " for i2 =  4:4+ nknots ;
end ;	%	end "if i < 3"

end ;	% 	end "for i:4+nknots"


end ; 	%	end "if smooth_spline_penalty == 1"

R=chol(xx + 10e-10*D) ;
B = inv(R') ;
[U,C] = eig(B*D*B') ;
Z = xm*B'*U ;

Zy = Z'*y ;
ZZ=Z'*Z ;

m = length(penwt) ; 
beta = zeros(size(xm,2),m) ;
%yhat = zeros(n,m) ;
asr = zeros(m,1) ;
gcv = asr ;
trsd = asr ;
trsdsd = asr ;
dfres = asr ;
ssy = y'*y ;

	for i=1:m ;
	oneld = 1 ./ (1 + penwt(i)*diag(C)) ;
	trsd(i) = sum(oneld ) ;
	alpha(:,i) = (Zy) .* oneld ;
	asr(i) =  (ssy - 2*Zy'*alpha(:,i) + alpha(:,i)'*ZZ*alpha(:,i))/n ;
	trsdsd(i) = sum(oneld.*oneld) ;
	dfres(i) = n -2*trsd(i) + trsdsd(i) ;
	gcv(i) = asr(i) / (1- gcvfact*trsd(i)/n)^2 ;
	end ;

imin = min(find(  (gcv==min(gcv)) ) ) ;
a = penwt(imin) ;

dffit= trsd ;
alpha = alpha(:,imin) ;
beta = B'* U * alpha ;
dbeta = length(beta) ;
yhat = Z*alpha ;
sigma2hat = n*asr(imin) ./ dfres(imin) ;

oneld = 1 ./ (1 + penwt(imin)*diag(C)) ;
postvaralpha = sigma2hat*diag(oneld);
postvaryhat = (Z.*(Z*postvaralpha))*ones(dbeta,1) ;
postvarbeta = B'*U*postvaralpha*U'*B ;

Z = [] ;

xmder = xm ;


xmder(:,1) = 0*xmder(:,1) ;

for i = 2:degree+1 ;
	xmder(:,i) = (i-1)*(abs(xm(:,i))).^( (i-2)/(i-1)  ) .* (sign(x)).^(i-2);
end ;

for i=degree+2:degree+1+nknots ;
	xmder(:,i) = degree* (abs(xm(:,i))).^(  (degree-1)/degree ) ...
		.* (xm(:,i) > 0);
end ;



yhatder = xmder*beta ;
postvaryhatder = ((xmder*postvarbeta).*xmder)*ones(dbeta,1) ;

clear xmder ;

if istd == 1 ;
	yhat = yhat + meany ;
end ;
ulimit = yhat + 2*sqrt(postvaryhat) ;
llimit = yhat - 2*sqrt(postvaryhat) ;
ulimitder = yhatder + 2*sqrt(postvaryhatder) ;
llimitder = yhatder - 2*sqrt(postvaryhatder) ;

xgrid = linspace(min(x),max(x),200)' ;
xgridm = powerbasis01(xgrid,degree,knots) ;

if istd == 1 ;
	xgrid = meanx + stdx*xgrid ;
	mhat = meany + xgridm*beta ;
	knots = meanx + stdx*knots ;
	yhdatder = yhatder./stdx ;
	ulimitder = ulimitder./stdx ;
	llimitder = llimitder./stdx ;

else ;
	mhat = xgridm*beta ;
end ;

if n < 5001 ;
fit = struct('yhat',yhat,'beta',beta,'gcv',gcv,'imin',imin,'dffit',dffit, ...
	'knots',knots,'postvarbeta',postvarbeta,'postvaryhat',postvaryhat, ...
	'xm',xm,'xx',xx,'a',a,'penwt',penwt,'sigma2hat',sigma2hat, ...
	'dfres',dfres, ...
	'yhatder',yhatder,'postvaryhatder',postvaryhatder, ...
	'ulimit',ulimit,'llimit',llimit,'ulimitder',ulimitder, ...
	'llimitder',llimitder, ...
	'asr',asr,'xgrid',xgrid,'mhat',mhat,'x',meanx+stdx*x,'degree',degree) ;
else ;
fit = struct('yhat',yhat,'beta',beta,'gcv',gcv,'imin',imin,'dffit',dffit, ...
	'knots',knots,'postvarbeta',postvarbeta,'postvaryhat',postvaryhat, ...
	'xx',xx,'a',a,'penwt',penwt,'sigma2hat',sigma2hat, ...
	'dfres',dfres, ...
	'yhatder',yhatder,'postvaryhatder',postvaryhatder, ...
	'ulimit',ulimit,'llimit',llimit,'ulimitder',ulimitder, ...
	'llimitder',llimitder, ...
	'asr',asr,'xgrid',xgrid,'mhat',mhat,'x',meanx+stdx*x,'degree',degree) ;

end ;
