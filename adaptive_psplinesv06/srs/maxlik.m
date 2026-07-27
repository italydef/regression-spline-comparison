function [betanew,varinfo,se] = maxlik(fun,beta,delta,x1,x2,x3,x4);
% Computes the MLE assuming independent data.
% Uses the Fisher method of scoring.
% fun is the loglikelihood at each observation.  fun must have
% argument (beta,x1,x2,x3,x4).  The x's are used to call in data.
%
% beta is the starting value.
% delta is the parameter for numerical derivatives.
% USE:  [betanew,varinfo,se] = maxlik(fun,beta,delta,x1,x2,x3,x4);
%
% Returns beta = MLE and varinfo = estimated cov. matrix using
%                and se = sqrt(diag(varinfo))
% an approximate Fisher information matrix.
%
%
m= size(beta,1) ;
fdiff = 1 ;
i =  1 ;
  while ( (i < 50) & (fdiff > .05) ) ;
  fbase = eval([fun,'(beta,x1,x2,x3,x4)']) ;
  n = size(fbase,1) ;
  fbase = sum(fbase) ; 
  gl = grad(fun,beta,delta,x1,x2,x3,x4) ; 
  h = inv(gl' * gl) ;
  del = h * (ones(1,n)*gl)' ;;
  betatest = beta + del ;
  fnew =  sum(eval([fun,'(betatest,x1,x2,x3,x4)']))  ;
    j = 1 ;
    while ( (fnew < fbase) & (j < 10) ) ;
    del = del/2 ;
    betatest = beta + del ;
    fnew =  sum(eval([fun,'(betatest,x1,x2,x3,x4)'])) ; 
    j = j + 1 ;
    end ;
  beta = betatest ;
  fdiff = fnew - fbase ; 
  i = i + 1 ;
  end ;
loglike = fnew ;
i 
loglike
fdiff
betanew = beta ;
varinfo = -inv(hessian(fun,beta,sqrt(delta),x1,x2,x3,x4)) ;
se =sqrt(diag(varinfo)) ;
'MLE     Std. Errors'
[betanew,se] 
