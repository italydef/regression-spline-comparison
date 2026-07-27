function [x,y,f] = rc_curve(n,j,errorvar,seed);

% This function generates data from the function given in Ruppert et al(1999)
% 
%   INPUTS
%
% n=number of data points
% 
% j=parameter determining spatial heterogeneity i.e. with
%     j=3 === small spatial heterogeneity
%     j=4 === Moderate spatial heterogeneity
%     j=6 === severe spatial heterogeneity
% 
% errorvar = adds the noise to the data generated
% 
% seed = seed of the random number generator
% 
%   OUPUTS
%   
%   x = x co-ordinate points linearly spaced over (0,1)
%   y = f + error
%   f = function evaluated at above x's 
%
%   Modified 9/24/02  
  
  randn('state',seed);  
  x = linspace(0, 1, n)';
  e = normrnd(0,sqrt(errorvar),n,1);
  numr=2*pi*(1+((2^((9-(4*j))/5))));
  denr=((2^((9-(4*j))/5)))+x;
  f = sqrt(x.*(1-x)).*(sin(numr./denr));
  y = f+e;
