function m = mmarron(x,j) ;
% mean function proposed by Steve Marron and used by Matt Wand in his study
% of regression splines
%
jterm = 2^((9-4*j)/5) ;
m = sqrt(x.*(1-x)) .* sin(2*pi*(1+jterm)./(x+jterm)) ;
