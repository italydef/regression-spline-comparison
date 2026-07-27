function [beta,stderror] = nlin2(fun,bstart,x,y,string) ;
%
%
%  Gauss-Newton nonlinear least squares algorithm                           
%  fun is the regression function and must be defined in an m-file somewhere               %  bstart = regression parameters starting values (m dimensional)
%  y = dependent variables  (n dimensional) 
%  x = independent variable vectors (n dimensional)                            
%  tol:  iteration stops if relative improvement falls below tol 
%        relative improvement is                                 
%        (ssold - ssnew)./max([ssnew tol])   
%  maxiter = maximum number of iterations  
%  iprint : set equal to 1 to print results of each iteration                  
%  CALLS these m-files: grad2
%
%  USAGE:   [beta,stderror] = nlin2(fun,bstart,x,y,string) ;
%
%  Last edited: October 10, 1997
%

tol = 1e-6 ;
maxiter = 25 ;
iprint = 0 ;
deltagrad=1e-6 ;   % delta for numerical gradients
%global iter sigma2 yhat residuals lambda;

eval('string') ;  %  Use string to modify the defaults

if size(bstart,1) == 1 ;
bstart=bstart';
end;

iter = 1 ;
crit = 2 * tol ;   % start crit > tol so we don't stop at the first step
beta=bstart ;

	while (crit > tol) & (iter <= maxiter) ;
	betaold = beta ;
	g =  grad2(fun,beta,x,deltagrad,string) ;  % compute gradient

	res = ( y - eval(fun) );
	ssold = sum( res.^2 ) ;
	delta = inv(g'*g)*g'*res ;  % delta is the Gauss-Newton updating
                                    % direction
	isqu = 1 ;              %   Begin squeezing   
	ok = 0 ;   % ok = 0 means we accepted the update
	    while (isqu <= 10) & (ok==0) ;
	    beta = betaold + delta ;
	    res = y - eval(fun) ;
	    ssnew = sum( (res).^2 ) ;
		if ssnew < ssold ;  %  if the criterion is met then update
		beta = betaold + delta ;     %  update
		ok = 1 ;            %  ok = 1 since we updated
		else ;
		delta = delta./2 ;  % half delta (called "squeezing")
                                    % if we don't update; then try
                                    % again
		isqu = isqu + 1 ;   % count the number of "squeeze"
		end;
	    end ;                  %   End squeezing   

	    if isqu == 11 ; crit = 0 ;  % Updata crit
	    else ;
	    crit = (ssold - ssnew) / max([tol; ssnew]) ; % we will stop if
                                     % crit is less than tol
	    end ;
	ssold = ssnew ;             % update the sum of squares
             if iprint == 1 ;
	    
	     '         iteration =  '; % Print summary of last iteration
	    iter 
	     '                 beta =  ';
	    beta'
	     'number of squeezes =  ';
	    isqu 
	     '    sum of squares = ';
	    ssnew
	     'stopping criterion =  ' ;
	    crit
            end ;
	iter = iter + 1 ;
	end ;
n = size(y,1) ;
sigma2 = sum(res.^2)./(n-size(beta,1)) ;
varmatrix = sigma2.* inv(g'*g) ;
stderror = sqrt(diag(varmatrix)) ;
yhat = eval(fun) ;
residuals = y - yhat ;
beta 
sigma2
stderror

