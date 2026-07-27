function y = hessian(fun,x,del,y1,y2,y3,y4) ;
%  Returns the hessian of fun
%  Usage: y = hess(fun,x,del,y1,y2,y3,y4) ;
%
n = rows(x) ;
ey = eye(n) ;
f1 = zeros(n,1) ;
f2 = zeros(n,n) ;
f = sum(eval([fun,'(x,y1,y2,y3,y4)']) );
  for i = 1:n ;
  f1(i) = sum(eval([fun,'(x+del*ey(:,i),y1,y2,y3,y4)'])) ;
    for j = 1:n ;
    f2(i,j) = sum(eval([fun,'(x+del*(ey(:,i)+ey(:,j)),y1,y2,y3,y4)'])) ;
      if j > i ; 
      f2(j,i) = f2(i,j) ; 
      end ;
    end ;
  end ;
y = (f2 - f1*ones(1,n) - ones(n,1)*f1' + f) / (del*del) ;
