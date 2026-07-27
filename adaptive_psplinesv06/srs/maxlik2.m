function [betanew,varinfo,se] = maxlik2(fun,beta,x,string);
%
%
%   Computes the MLE assuming independent data.
%   Uses the Fisher method of scoring.
%   fun evaluates the loglikelihood at each observation (string)
%   delta is the parameter for numerical derivatives.
%   USASE:  [betanew,varinfo,se] = maxlik(fun,beta,x,string);
%
%   Returns beta = MLE and varinfo = estimated cov. matrix using
%                  an approximate Fisher information matrix.
%                  and se = sqrt(diag(varinfo))
%
%   CALLS:  grad2, hessian2
%
delta = 10e-5 ;    %  Set defaults
eval(string) ;     % Can be used to change defaults
m= rows(beta) ;
fdiff = 1 ;        %  fdiff is change in log-likelihood
i =  1 ;           %  i is iteration number

  while ( (i < 50) & (fdiff > .01) ) ;  
  fbase = eval(fun) ;     %  Get current log-likelhood
  n = rows(fbase) ;
  fbase = sum(fbase) ; 
  gl = grad2(fun,beta,x,delta) ; 
  h = inv(gl' * gl) ;
  del = h * (ones(1,n)*gl)' ;    %  Direction to increase log-likelihood
  betaold = beta ;
  beta = betaold + del ;         %  New MLE
  fnew =  sum(eval(fun))  ;
    j = 1 ;
    while ( (fnew < fbase) & (j < 10) ) ;   %  Check is new MLE improves
    del = del/2 ;                %  Half step if new MLE doesn't improve        
    beta = betaold + del ;
    fnew =  sum(eval(fun)) ; 
    j = j + 1 ;
    end ;
  fdiff = fnew - fbase ;      %  Calculate increase in log-likelihood
  i = i + 1 ;
  end ;
loglike = fnew ;
%  print out results
['number of iterations = ',num2str(i)]
loglike
fdiff
varinfo = -inv(hessian2(fun,beta,x,sqrt(delta))) ;  % inverse Fisher info
se =sqrt(diag(varinfo)) ;
'MLE     Std. Errors'
[beta,se] 
betanew = beta ;
