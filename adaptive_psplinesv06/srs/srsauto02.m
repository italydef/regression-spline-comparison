function fit = srsauto02(x,y,degree,penwt,maxnknots,boundstab, ...
	smooth_spline_penalty,gcvfact) ;

%
%	Calls PsplineDR01 with successively larger values of nknots.
%	Stops when gcv is not reduced by at
%	least 2%.
%	
%	Values of nknots tried:  5, 10, 20, 40, 80, 120
%
%
%		INPUT (required)
%	x = independent variable (n by 1)
%	y = dependent variable (n by 1)
%
%		INPUT (optional)
%	degree = spline degree (default = 2)
%	penwt = vector of possible penalty weights 
%		(default = logspace(-10,10,51)
%	maxnknots = maximum number of knots 
%	boundstab =  boundary stability parameter (default = 0)
%
% 		RETURNS
%	Returns a structure with the following fields -
%
%	fit1 = fit with chosen knots
%	nknotschosen = number of knots chosen
%
%
%
%	Last edit: 8/16/2001
%
%	Calls: PsplineDR03

if (nargin < 8 | isempty(gcvfact) == 1 ) ;
	gcvfact = 1 ;
end ;

if (nargin < 7 | isempty(smooth_spline_penalty) == 1 ) ;
	smooth_spline_penalty = 0 ;
end ;

if (nargin < 6 | isempty(boundstab) == 1 );
	boundstab = 0 ;
end ;

if (nargin < 5 | isempty(maxnknots) == 1) ;
	maxnknots = 120 ;
end ;

if (nargin < 4 | isempty(penwt) == 1) ;
	penwt = logspace(-10,10,51) ;
end ;

if (nargin < 3 | isempty(degree) == 1) ;
	degree = 2 ;
end ;


n = length(x) ;
nu = length(unique(x)) ;

nknots = [5 10 20 40 80 120]' ;
nknots = nknots( (nknots < min([maxnknots+1  nu - degree - 1]) ) ) ;


fit1 = [] ;
fit2 = [] ;
fit3 = [] ;
fit4 = [] ;
fit5 = [] ;
fit6 = [] ;


for i=1:length(nknots) ;

fit = PsplineDR03(x,y,struct('degree',degree,'nknots',nknots(i), ...
	'smooth_spline_penalty',smooth_spline_penalty,'gcvfact',gcvfact) ) ;

if i == 1 ;
fit1 = fit ;
elseif i == 2 ;
fit2 = fit ;
elseif i == 3 ;
fit3 = fit ;
elseif i == 4 ;
fit4 = fit ;
elseif i == 5 ;
fit5 = fit ;
elseif i == 6 ;
fit6 = fit ;
end ;

gcv(i) = fit.gcv(fit.imin) ;
end ;

iminknots = min( find(min(gcv)==gcv)) ;

if iminknots == 1 ;
fitauto = fit1 ;
elseif iminknots == 2 ;
fitauto = fit2 ;
elseif iminknots == 3 ;
fitauto = fit3 ;
elseif iminknots == 4 ;
fitauto = fit4;
elseif iminknots == 5 ;
fitauto = fit5 ;
elseif iminknots == 6 ;
fitauto = fit6 ;
end ;

nknotschosen = nknots(iminknots) ;

fit = struct('fit1',fit1,'fit2',fit2,'fit3',fit3,'fit4',fit4, ...
	'fit5',fit5,'fit6',fit6, ...
	'fitauto',fitauto,'iminknots',iminknots,'nknotschosen',nknotschosen) ;
