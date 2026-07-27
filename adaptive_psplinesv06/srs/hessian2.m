function y = hessian2(fun,beta,x,del) ;
%  Returns the hessian of fun
%  Usage: y = hessian2(fun,beta,x,del) ;
%
m = rows(beta) ;
ey = eye(m) ;
f1 = zeros(m,1) ;
f2 = zeros(m,m) ;
f = sum(eval(fun) );
betaold = beta ;
  for i = 1:m ;
  f1(i) = sum(eval(fun)) ;
    for j = 1:m ;
    beta = betaold + del*(ey(:,i)+ey(:,j)) ;
    f2(i,j) = sum(eval(fun)) ;
      if j > i ; 
      f2(j,i) = f2(i,j) ; 
      end ;
    end ;
  end ;
y = (f2 - f1*ones(1,m) - ones(m,1)*f1' + f) / (del*del) ;
