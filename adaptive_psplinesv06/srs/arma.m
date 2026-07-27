function [arest,se,t,e,th]=arma(y,ndiff,nar,nma,center);
%	Computes arma model estimates by calling armax of the 
%	identifcation toolbox
%	
%		INPUT (required)
%	y	= time series (univariate)
%
%		INPUT (optional)
%	ndiff	= number of times the series should be differenced (default
%			= 0)
%	center	= (=1 if the series should be centered) (default = 1)
if nargin < 2;
ndiff = 0 ;
end ;

if nargin < 3 ;
center = 1 ;
end ;

y=y-mean(y) ;
npar = nar + nma ;
th = armax(y,[nar nma]) ;
cov = th(4:3+npar,1:npar) ;
se=sqrt(diag(cov))' ;
arest = th(3,1:npar) ;
t = arest./se ;
disp(['std errors =       ',num2str(se)]);
disp(['ARMA estimates =   ', num2str(arest)] )
disp(['t statistics =     ',num2str(t)]) ;
e=resid(y,th) ;
