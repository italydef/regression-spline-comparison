function gamvars = randgam(nrow, ncol, a, b)
%
% randgam(nrow, ncol, a, b)
%    returns a matrix of dimension nrow by ncol of
%    observations from a gamma distribution with parameters
%    a and b.  (i.e. density is proportional to  x^(a-1) exp(-b*x)
%
%    Currently works only if a is an integer.
%
%    Written by Elizabeth Slate
%
N = nrow * ncol;
U = rand(a, N);
E = -log(U)/b;
G = sum(E); 
gamvars = reshape(G,nrow,ncol);
end



