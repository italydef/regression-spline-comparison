function y = pdfn(x) ;
% Returns probability density function of the N(0,1) distribution
% USE: y = pdfn(x)
y = exp(- .5 * x.^2) / (sqrt(2*pi)) ; 
