% 
% Code generating the fits
% Modified 6/02/02
% 
clear all;
n=400; % number of data points
errorvar=0.04 ; % error variance
seed = 0 ; % seed for random number 
nIter =100000; % number of Iterations for MCMC
burn_in = 50000; % burn-in for MCMC
degree=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% j=6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




[x,y,f] = rc_curve(n,6,errorvar,seed);


nknots=[40] ; % number of knots
nsubknots=[4];% number of subknots



   
    fit=bayes_pspline(x,y,nknots,nsubknots,degree,nIter,burn_in) ;
    mse = ( mean((f-fit.yhat).^2));
    fit.mse = mse;
    
      
    
 display([' ******  '])   
 display([' ******  Estimated MSE =  ', num2str(mse)])
 display([' ******  ']) 
        
        
   