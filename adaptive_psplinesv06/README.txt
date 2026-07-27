

SPATIALLY ADAPTIVE BAYESIAN PENALIZED REGRESSION SPLINES (P-SPLINES) ALGORITHM (v6.0)


from paper by Baladandayuthapani, Mallick and Carroll (JCGS) (2005)

This is the algorithm used in the paper. Acknowledgements to


Veerabhadran Baladandayuthapani, Bani K. Mallick and Raymond J. Carroll,
Department of Statistics, Texas A&M University, College Station, Texas, USA.


Any problems/comments email: veera@stat.tamu.edu



*************************************************************************

CONTENTS:

1. bayes_pspline.m 		    : for univariate functions
2. example.m 	   		    : example call for univariate functions
3. srs directory	            : David Ruppert's Matlab codes from
				http://www.orie.cornell.edu/~davidr/				    
	

*************************************************************************

INSTRUCTIONS:

Assuming that the programs have already been unpacked. The codes use starting values 
generated from the algorithm given in Ruppert and Carroll (2000). The MATLAB codes 
can be downloaded from  http://www.orie.cornell.edu/~davidr/  

 For Univariate functions.

	The function bayes_psline() fits a smoothed regression spline using a P-spline 
	basis. 
	
	For example SEE : example.m


To replicate results in the Program:

Use knots = specify number of knots
    nsubknots = specify number of subknots
    nIter = number of MCMC iterations
    burin_in = burn-in for MCMC
    



ABOUT THE PROGRAM

1. The number of iterations and the length of burn-in used vary widely from problem
   to problem, hence various criteria for convergence has been suggested, but has not been
   implemented here. In general after 10,000 burn-in iterations the sampler would have
   almost converged but most problems may require less iterations.


2. We have not dealt into knot selection here but it is imperative that a certain minimum number 
   of knots are chosen to capture the spatial variability of the function. Too few knots may 
   smooth out important peaks and too many knots may make the fit more spiky.


(You are welcome to email be with any nice modifications of the program)

*************************************************************************

COPYRIGHT

(c) Copyright Veerabhadran Baladandayuthapani (2006)

You are free to use this program, for non-commercial purposes only 
under two conditions:

1. This note is note removed
2. Publications using this code should cite "Spatially Adaptive Bayesian Penalized Regression 
   Splines (P-splines)" paper.


*************************************************************************

modified 4/4/06




	





