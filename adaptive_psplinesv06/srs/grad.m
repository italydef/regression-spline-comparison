function g = grad(fun,theta,deltagrad,x1,x2,x3,x4) ;
% Computes the grad of 'fun'
% Usage:  g = grad('fun',theta,delta,x1,x2,x3,x4,x5) ;
%
m = rows(theta) ;
if m== 1 ; theta = theta' ; end ;
m=rows(theta) ;

ii = eye(m) ;
fbase = eval([fun,'(theta,x1,x2,x3,x4)']) ;

for i = 1:m ;
theta2 = theta + deltagrad * ii(:,i) ;
g(:,i) = (eval([fun,'(theta2,x1,x2,x3,x4)']) - fbase) / deltagrad ;
end ;

