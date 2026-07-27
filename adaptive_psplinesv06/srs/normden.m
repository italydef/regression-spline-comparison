function [d,logd] = normden(x,mu,cov,icov) ;
%
%	Computes the multivariate normal density with mean = mu
%	and covariance matrix = cov and icov = inverse of cov
%
%	icov is optional
%
%	Last edit: 11/23/97
%
%
if nargin < 4 ;
icov = inv(cov) ;
end ;

L = length(x) ;
quad = (x-mu)'*icov*(x-mu) ;
expo = exp(-.5*quad) ;
d = (2*pi)^(-L/2) * det(cov)^(-1/2) * expo;
logd = -L/2*log(2*pi) + .5*log(det(cov)) -.5*quad ;
