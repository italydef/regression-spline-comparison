% 	Program to test loccv 
%
n = 200 ;
power = 2 ;
ism = 1 ;
if ism == 1 ;
	x=rand(n,1) ;
	x = sort(x) ;
	y= sin(15*x.^power) + randn(n,1)/5 ;
tic
	[yhatloccv,spansloccv,hloccv,yhatcv,spancv,hcv] = ...
		loccv(x,y,2) ;
toc
end ;
[yhatgl,yhatloc] = srslocal(x,y,2) ;

subplot(2,1,1) ;
plot(x,yhatcv,'b-.',x,yhatloccv,'r-',x,y,'*',x,yhatloc,'g--') ;
legend('CV', 'loc CV', 'data','local srs') ;

subplot(2,1,2) ;
plot(x,spancv*ones(n,1),'b-.',x,spansloccv,'r-') ;
legend('CV span', 'local CV span') ;
