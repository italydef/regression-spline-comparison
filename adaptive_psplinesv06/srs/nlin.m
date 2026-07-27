function beta = nlin(fun,bstart,x1,x2,x3,x4,y,tol,maxiter,iprint) ;


%  Gauss-Newton nonlinear least squares algorithm                           
%  fun is the regression function and must be defined in an m-file somewhere               %  bstart = regression parameters starting values (m dimensional)
%  y = dependent variables  (n dimensional) 
%  x1,x2,x3,x4 = four independent variable vectors (each in n dimensional)                                    %  If only some of the x's are needed, the others can be set to 0
%  tol:  iteration stops if relative improvement falls below tol 
%        relative improvement is                                 
%        (ssold - ssnew)./maxc(ssnew|tol)   
%  maxiter = maximum number of iterations  
%  iprint : set equal to 1 to print results of each iteration                  
%  CALLS these m-files: rows, grad
%

if rows(bstart) == 1 ;
bstart=bstart';
end;

deltagrad=1e-6 ;   % delta for numerical gradients
iter = 1 ;
crit = 2 * tol ;   % start crit > tol so we don't stop at the first step
b=bstart ;

	while (crit > tol) & (iter <= maxiter) ;
	g =  grad(fun,b,deltagrad,x1,x2,x3,x4) ;  % compute gradient
	res = ( y - eval([fun,'(b,x1,x2,x3,x4)']) );
	ssold = sum( res.^2 ) ;
	delta = inv(g'*g)*g'*res ;  % delta is the Gauss-Newton updating
                                    % direction
	    isqu = 1 ;              %   Begin squeezing   
	    ok = 0 ;   % ok = 0 means we accepted the update
	    while (isqu <= 10) & (ok==0) ;
	    res = (y - eval(  [fun,'((b + delta),x1,x2,x3,x4)'])  ) ;
	    ssnew = sum( (res).^2 ) ;
		if ssnew < ssold ;  %  if the criterion is met then update
		b = b + delta ;     %  update
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
	     '                 b =  ';
	    b'
	     'number of squeezes =  ';
	    isqu 
	     '    sum of squares = ';
	    ssnew
	     'stopping criterion =  ' ;
	    crit
            end ;
	iter = iter + 1 ;
	end ;

beta = b ;  % Final value of b that is output
